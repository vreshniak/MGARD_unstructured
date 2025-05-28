/*
 * make_ordering.cpp :
 *
 *  Created on: October 2, 2024
 *      Author: Viktor Reshniak
 *
 *	Reorder serialized DOFs
 */


#include <adios2.h>
#include <vector>
#include <algorithm>    // std::find

#include <iostream>
#include <fstream>

#include <mpi.h>
#include <chrono>

#include "serialization.hpp"

#include <nlohmann/json.hpp>
using json = nlohmann::json;


bool verbose = true;


// /**
//  * Calculate global node connectivity
//  */
// template <typename Container>
// auto ComputeGlobalConnectivity( const Container ElementConnectivity, size_t num_nodes, int nodes_per_el )
// -> std::vector<std::set<size_t>>
// {
// 	std::vector<std::set<size_t>> GlobalConnectivity(num_nodes);

// 	// loop through elements
// 	for(size_t k = 0; k < ElementConnectivity.size(); k+=nodes_per_el){
// 		// loop through nodes in each element
// 		for(int i = 0; i < nodes_per_el; i++){
// 			// index of current element
// 			size_t n1 = ElementConnectivity[k+i];
// 			// another loop through nodes in same element
// 			for(size_t j = 0; j < nodes_per_el; j++){
// 				// index of current element
// 				size_t n2 = ElementConnectivity[k+j];
// 				// update connected neighbors of each node in the element
// 				if (n1!=n2){
// 					GlobalConnectivity[n1].insert(n2);
// 					GlobalConnectivity[n2].insert(n1);
// 				}
// 			}
// 		}
// 	}
// 	return GlobalConnectivity;
// }


// template <typename T>
// auto reorder(const std::vector<T> &vec, const std::vector<size_t> &order)
// -> std::vector<T>
// {
// 	std::vector<T> new_vec(vec.size());

// 	for (size_t i=0; i<vec.size(); i++)
// 		new_vec[i] = vec[order[i]];

// 	return new_vec;
// }


// template <typename Container>
// auto find_node_order(const std::vector<std::set<size_t>> &GlobalConnectivity, const Container &coos)
// -> std::vector<int64_t>
// {
// 	size_t num_nodes = GlobalConnectivity.size();

// 	std::vector<bool>   node_notvisited(num_nodes, true);
// 	std::vector<int64_t> node_order(num_nodes, 0);
// 	node_notvisited[0] = false;
// 	size_t curr_node = 0;
// 	size_t node_num_visited = 1;
// 	while (node_num_visited<num_nodes){
// 		// find connected nodes that have not been visited yet
// 		std::vector<size_t> nns;
// 		for (const auto & n : GlobalConnectivity[curr_node])
// 			if (node_notvisited[n])
// 				nns.push_back(n);

// 		// if((node_num_visited-1)%1000==0){
// 		// 	std::cout << curr_node << ": ";
// 		// 	for (const auto & n : nns) std::cout << n << " ";
// 		// 	std::cout << std::endl;
// 		// }

// 		size_t next_node;
// 		if (nns.size()>0){
// 			double min_dst = 1.e15;
// 			for (size_t i=0; i<nns.size(); i++){
// 				size_t candidate = nns[i];
// 				double dst = 0;
// 				for (auto &coo : coos)
// 					dst += (coo[curr_node]-coo[candidate]) * (coo[curr_node]-coo[candidate]);
// 				if (dst<min_dst){
// 					min_dst   = dst;
// 					next_node = candidate;
// 				}
// 			}
// 		}else{
// 			// find first not visited node
// 			next_node = std::find(node_notvisited.begin(), node_notvisited.end(), true) - node_notvisited.begin();
// 		}
// 		node_notvisited[next_node]   = false;
// 		node_order[node_num_visited] = next_node;

// 		curr_node = next_node;
// 		node_num_visited++;

// 		// if(node_num_visited%1000==0)
// 		// 	std::cout << node_num_visited << " " << num_nodes << " " << curr_node << std::endl;
// 	}

// 	return node_order;
// }



int main(int argc, char *argv[])
{
	///////////////////////////////////////////////////////////////////////////
	// MPI

	int rank, nranks;
#if ADIOS2_USE_MPI
	int provided;

	// MPI_THREAD_MULTIPLE is only required if you enable the SST MPI_DP
	MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nranks);
#else
	rank   = 0;
	nranks = 1;
#endif


	///////////////////////////////////////////////////////////////////////////
	// config

	// read config file
	if (argc<4){
		std::cerr << "Usage: make_ordering config.json input.bp output.bp, `make_ordering` got " << argc-1 << " parameters" << std::endl;
		return -1;
	}
	std::string config_file_name = argv[1];
	std::string input_file_name  = argv[2];
	std::string output_file_name = argv[3];

	std::ifstream f(config_file_name);
	json config = json::parse(f);

	// variable names
	std::vector<std::string> coo_names = config["coo_names"];
	std::string connectivity_name = config["connectivity"];
	std::string node_order_name = config["node_order"];


	// setup ADIOS readers/writers
#if ADIOS2_USE_MPI
	adios2::ADIOS adios(MPI_COMM_WORLD);
#else
	adios2::ADIOS adios;
#endif
	adios2::IO reader_io = adios.DeclareIO("BPReader");
	adios2::IO writer_io = adios.DeclareIO("BPWriter");
	reader_io.SetEngine("BPFile");
	writer_io.SetEngine("BPFile");
	adios2::Engine bpReader = reader_io.Open(input_file_name,  adios2::Mode::Read);
	adios2::Engine bpWriter = writer_io.Open(output_file_name, adios2::Mode::Write);


	///////////////////////////////////////////////////////////////////////////
	// define output ADIOS variables

	// connectivity
	adios2::Variable<int64_t> out_adios_connectivity = writer_io.DefineVariable<int64_t>(connectivity_name, {}, {}, {adios2::UnknownDim});

	// coordinates
	std::vector<adios2::Variable<double>> out_adios_coos(coo_names.size());
	for (int i=0; i<coo_names.size(); i++)
		out_adios_coos[i] = writer_io.DefineVariable<double>(coo_names[i], {}, {}, {adios2::UnknownDim});

	// node ordering
	adios2::Variable<int64_t> out_adios_order = writer_io.DefineVariable<int64_t>(node_order_name, {}, {}, {adios2::UnknownDim});


	///////////////////////////////////////////////////////////////////////////


	double total_time = 0;


	// time stepping
	while (true){

		adios2::StepStatus read_status = bpReader.BeginStep(adios2::StepMode::Read, -1.0f);
		if (read_status!=adios2::StepStatus::OK) break;
		size_t step = bpReader.CurrentStep();
		bpWriter.BeginStep();


		///////////////////////////////////////////////////////////////////////////
		// define input ADIOS variables

		adios2::Variable<int64_t> adios_connectivity = reader_io.InquireVariable<int64_t>(connectivity_name);

		std::vector<adios2::Variable<double>> adios_coos(coo_names.size());
		for (int i=0; i<coo_names.size(); i++)
			adios_coos[i] = reader_io.InquireVariable<double>(coo_names[i]);

		///////////////////////////////////////////////////////////////////////////
		// find reordering

		// number of blocks/subdomains
		auto blocks = bpReader.BlocksInfo(adios_connectivity, 0);
		size_t nblocks = blocks.size();

		int block_id = 0;
		for (auto &block : blocks){
			auto start = std::chrono::high_resolution_clock::now();
			if (verbose) std::cout << "step = " << step << ", blockID = " << block.BlockID << " out of " << nblocks << std::endl;

			// allocate memory for reading variables
			std::vector<int64_t> ElementConnectivity;
			std::vector<std::vector<double>> coos(coo_names.size());

			// set block/subregion for the variables
			adios_connectivity.SetBlockSelection(block.BlockID);
			for (auto &coo : adios_coos) coo.SetBlockSelection(block.BlockID);

			// read variables
			bpReader.Get(adios_connectivity, ElementConnectivity, adios2::Mode::Sync);
			for (int i=0; i<coo_names.size(); i++) bpReader.Get(adios_coos[i], coos[i], adios2::Mode::Sync);
			bpReader.PerformGets();


			//////////////////////////////////////////////////////////////////////////////////
			// reordering

			size_t num_nodes = coos[0].size();

			// compute connectivity of the node graph
			auto GlobalConnectivity = ComputeGlobalConnectivity(ElementConnectivity, num_nodes, config["n_nodes_in_element"]);

			auto node_order = find_node_order(GlobalConnectivity, coos);

			// node_order  = [4, 3, 0, 2, 1]
			// reorder_map = [2, 4, 3, 1, 0], i.e.,
			// 0 -> 2
			// 1 -> 4
			// 2 -> 3
			// 3 -> 1
			// 4 -> 0
			// std::vector<size_t> reorder_map(node_order.size());
			// for (size_t i=0; i<reorder_map.size(); i++)
			// 	reorder_map[node_order[i]] = i;

			// for (size_t i=0; i<ElementConnectivity.size(); i++)
			// 	ElementConnectivity[i] = reorder_map[ElementConnectivity[i]];


			//////////////////////////////////////////////////////////////////////////////////
			// save ordering

			// connectivity
			out_adios_connectivity.SetSelection(adios2::Box<adios2::Dims>({}, {ElementConnectivity.size()}));
			bpWriter.Put<int64_t>(out_adios_connectivity, ElementConnectivity.data(), adios2::Mode::Sync);

			// coordinates
			for (int i=0; i<coo_names.size(); i++){
				out_adios_coos[i].SetSelection(adios2::Box<adios2::Dims>({}, {coos[i].size()}));
				bpWriter.Put<double>(out_adios_coos[i], coos[i].data(), adios2::Mode::Sync);
			}

			// node ordering
			out_adios_order.SetSelection(adios2::Box<adios2::Dims>({}, {node_order.size()}));
			bpWriter.Put<int64_t>(out_adios_order, node_order.data(), adios2::Mode::Sync);


			bpWriter.PerformPuts();
			///////////////////////////////////////////////////////////////

			auto end = std::chrono::high_resolution_clock::now();
			auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
			total_time += (double)duration.count() / 1e6;
			if (verbose) std::cout << "\ttime = " << (double)duration.count() / 1e6 << " sec" << std::endl;


			///////////////////////////////////////////////////////////////

			if (config["max_block_id"]>=0 and block_id++>=config["max_block_id"]) break;
		}

		bpReader.EndStep();	// end logical step
		bpWriter.EndStep();	// end logical step

		if (step>=config["max_step"]) break;
	}

	if (verbose) std::cout << "Total time = " << std::setprecision(3) << total_time << std::endl;

	bpReader.Close();	// close engine
	bpWriter.Close();	// close engine

#if ADIOS2_USE_MPI
	MPI_Finalize();
#endif
	return 0;
}
