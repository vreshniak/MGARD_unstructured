#include <adios2.h>
#include <vector>
#include <algorithm>    // std::find

#include <iostream>
#include <fstream>

#include <mpi.h>
#include <chrono>

#include <nlohmann/json.hpp>

using json = nlohmann::json;


/**
 * Calculate global node connectivity
 */
template <typename Container>
auto ComputeGlobalConnectivity( const Container ElementConnectivity, size_t num_nodes, int nodes_per_el )
-> std::vector<std::set<size_t>>
{
	std::vector<std::set<size_t>> GlobalConnectivity(num_nodes);

	for(size_t k = 0; k < ElementConnectivity.size(); k+=nodes_per_el){
		for(int i = 0; i < nodes_per_el; i++){
			size_t n1 = ElementConnectivity[k+i];
			for(size_t j = 0; j < nodes_per_el; j++){
				size_t n2 = ElementConnectivity[k+j];
				if (n1!=n2){
					GlobalConnectivity[n1].insert(n2);
					GlobalConnectivity[n2].insert(n1);
				}
			}
		}
	}
	return GlobalConnectivity;
}


template <typename T>
auto reorder(const std::vector<T> &vec, const std::vector<size_t> &order)
-> std::vector<T>
{
	std::vector<T> new_vec(vec.size());

	for (size_t i=0; i<vec.size(); i++)
		new_vec[i] = vec[order[i]];

	return new_vec;
}


template <typename Container>
auto find_node_order(const std::vector<std::set<size_t>> &GlobalConnectivity, const Container &coos)
-> std::vector<size_t>
{
	size_t num_nodes = GlobalConnectivity.size();

	std::vector<bool>   node_notvisited(num_nodes, true);
	std::vector<size_t> node_order(num_nodes, 0);
	node_notvisited[0] = false;
	size_t curr_node = 0;
	size_t node_num_visited = 1;
	while (node_num_visited<num_nodes){
		// find connected nodes that have not been visited yet
		std::vector<size_t> nns;
		for (const auto & n : GlobalConnectivity[curr_node])
			if (node_notvisited[n])
				nns.push_back(n);

		// if((node_num_visited-1)%1000==0){
		// 	std::cout << curr_node << ": ";
		// 	for (const auto & n : nns) std::cout << n << " ";
		// 	std::cout << std::endl;
		// }

		size_t next_node;
		if (nns.size()>0){
			double min_dst = 1.e15;
			for (size_t i=0; i<nns.size(); i++){
				size_t candidate = nns[i];
				double dst = 0;
				for (auto &coo : coos)
					dst += (coo[curr_node]-coo[candidate]) * (coo[curr_node]-coo[candidate]);
				if (dst<min_dst){
					min_dst   = dst;
					next_node = candidate;
				}
			}
		}else{
			// find first not visited node
			next_node = std::find(node_notvisited.begin(), node_notvisited.end(), true) - node_notvisited.begin();
		}
		node_notvisited[next_node]   = false;
		node_order[node_num_visited] = next_node;

		curr_node = next_node;
		node_num_visited++;

		// if(node_num_visited%1000==0)
		// 	std::cout << node_num_visited << " " << num_nodes << " " << curr_node << std::endl;
	}

	return node_order;
}



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
	rank = 0;
	nranks = 1;
#endif


	///////////////////////////////////////////////////////////////////////////
	// config

	// read config file
	if (argc<4){
		std::cerr << "Config JSON file is required as input" << std::endl;
		return -1;
	}
	std::ifstream f(argv[1]);
	json config = json::parse(f);


	// variable names
	std::vector<std::string> coo_names = config["coo_names"];
	std::vector<std::string> var_names = config["var_names"];
	std::string connectivity_name = config["connectivity"];


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
	adios2::Engine bpReader = reader_io.Open(argv[2], adios2::Mode::Read);
	adios2::Engine bpWriter = writer_io.Open(argv[3], adios2::Mode::Write);


	// // ADIOS compression
	// adios2::Operator op = adios.DefineOperator("mgardplus", "mgardplus");


	///////////////////////////////////////////////////////////////////////////
	// define output ADIOS variables
	// adios2::Variable<size_t> inv_map_out = writer_io.DefineVariable<size_t>("inv_map", {}, {}, {adios2::UnknownDim});

	adios2::Variable<int64_t> out_adios_connectivity = writer_io.DefineVariable<int64_t>(connectivity_name, {}, {}, {adios2::UnknownDim});

	std::vector<adios2::Variable<double>> out_adios_coos(coo_names.size());
	for (int i=0; i<coo_names.size(); i++)
		out_adios_coos[i] = writer_io.DefineVariable<double>(coo_names[i], {}, {}, {adios2::UnknownDim});

	std::vector<adios2::Variable<double>> out_adios_vars(var_names.size());
	for (int i=0; i<var_names.size(); i++)
		out_adios_vars[i] = writer_io.DefineVariable<double>(var_names[i], {}, {}, {adios2::UnknownDim});

	///////////////////////////////////////////////////////////////////////////


	double total_time = 0;


	// time stepping
	while (true){


	adios2::StepStatus read_status = bpReader.BeginStep(adios2::StepMode::Read, -1.0f);
	size_t step = bpReader.CurrentStep();
	if (step>config["max_step"] or read_status!=adios2::StepStatus::OK)
		break;
	bpWriter.BeginStep();


	///////////////////////////////////////////////////////////////////////////
	// define input ADIOS variables

	adios2::Variable<int64_t> adios_connectivity = reader_io.InquireVariable<int64_t>(connectivity_name);

	std::vector<adios2::Variable<double>> adios_coos(coo_names.size());
	for (int i=0; i<coo_names.size(); i++)
		adios_coos[i] = reader_io.InquireVariable<double>(coo_names[i]);

	std::vector<adios2::Variable<double>> adios_vars(var_names.size());
	for (int i=0; i<var_names.size(); i++)
		adios_vars[i] = reader_io.InquireVariable<double>(var_names[i]);


	///////////////////////////////////////////////////////////////////////////
	// reorder


	// number of blocks/subdomains
	auto blocks = bpReader.BlocksInfo(adios_connectivity, 0);
	size_t nblocks = blocks.size();

	int block_id = 0;
	for (auto &block : blocks){
		auto start = std::chrono::high_resolution_clock::now();
		std::cout << "step = " << step << ", blockID = " << block.BlockID << " out of " << nblocks << std::endl;

		// allocate memory for reading variables
		std::vector<int64_t> ElementConnectivity;
		std::vector<std::vector<double>> coos(coo_names.size());
		std::vector<std::vector<double>> vars(var_names.size());

		// set block/subregion for the variables
		adios_connectivity.SetBlockSelection(block.BlockID);
		for (auto &coo : adios_coos) coo.SetBlockSelection(block.BlockID);
		for (auto &var : adios_vars) var.SetBlockSelection(block.BlockID);

		// read variables
		bpReader.Get(adios_connectivity, ElementConnectivity, adios2::Mode::Sync);
		for (int i=0; i<coo_names.size(); i++) bpReader.Get(adios_coos[i], coos[i], adios2::Mode::Sync);
		for (int i=0; i<var_names.size(); i++) bpReader.Get(adios_vars[i], vars[i], adios2::Mode::Sync);
		bpReader.PerformGets();


		//////////////////////////////////////////////////////////////////////////////////
		// reordering

		size_t num_nodes = coos[0].size();

		// compute connectivity of the node graph
		auto GlobalConnectivity = ComputeGlobalConnectivity(ElementConnectivity, num_nodes, config["n_nodes_in_element"]);

		// // print GlobalConnectivity
		// for (int curr_node=0; curr_node<GlobalConnectivity.size(); curr_node++){
		// 	std::cout << curr_node << " ";
		// 	for (const auto & n : GlobalConnectivity[curr_node])
		// 		std::cout << n << " ";
		// 	std::cout << std::endl;
		// }
		// return 0;

		auto node_order = find_node_order(GlobalConnectivity, coos);

		// for (int i=0; i<10; i++)
		// 	std::cout << node_order[i] << " ";
		// std::cout << std::endl;

		std::vector<size_t> reorder_map(node_order.size());
		for (size_t i=0; i<reorder_map.size(); i++)
			reorder_map[node_order[i]] = i;

		for (size_t i=0; i<ElementConnectivity.size(); i++)
			ElementConnectivity[i] = reorder_map[ElementConnectivity[i]];


		//////////////////////////////////////////////////////////////////////////////////
		// save reordered variables

		// double abs_tol = 1.e-3;
		// RHE_out.AddOperation(op, {{"accuracy", std::to_string(abs_tol)}, {"mode", "ABS"}});

		// set block/subregion for the variables
		out_adios_connectivity.SetSelection(adios2::Box<adios2::Dims>({}, {ElementConnectivity.size()}));
		bpWriter.Put<int64_t>(out_adios_connectivity, ElementConnectivity.data(), adios2::Mode::Sync);
		for (int i=0; i<coo_names.size(); i++){
			out_adios_coos[i].SetSelection(adios2::Box<adios2::Dims>({}, {coos[i].size()}));
			bpWriter.Put<double>(out_adios_coos[i], reorder(coos[i],node_order).data(), adios2::Mode::Sync);
		}
		for (int i=0; i<var_names.size(); i++){
			out_adios_vars[i].SetSelection(adios2::Box<adios2::Dims>({}, {vars[i].size()}));
			bpWriter.Put<double>(out_adios_vars[i], reorder(vars[i],node_order).data(), adios2::Mode::Sync);
		}

		// inv_map_out.SetSelection(adios2::Box<adios2::Dims>({}, {serial_inv_map.size()})); bpWriter.Put<size_t>(inv_map_out, serial_inv_map.data(), adios2::Mode::Sync);
		bpWriter.PerformPuts();


		///////////////////////////////////////////////////////////////

		auto end = std::chrono::high_resolution_clock::now();
		auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

		total_time += (double)duration.count() / 1e6;

		std::cout << "\ttime = " << (double)duration.count() / 1e6 << " sec" << std::endl;

		///////////////////////////////////////////////////////////////
		// print original and reordered variables

		// for (size_t i=0; i<50000; i++)
		// 	std::cout << vars[0][i] << " ";
		// std::cout << std::endl;

		// std::vector<double> reordered_var = reorder(vars[0],node_order);
		// for (size_t i=0; i<50000; i++)
		// 	std::cout << reordered_var[i] << " ";
		// std::cout << std::endl;

		// return 0;

		///////////////////////////////////////////////////////////////
		// print original and reordered coordinates

		// for (int coo=0; coo<3; coo++){
		// 	for (size_t i=0; i<50000; i++)
		// 		std::cout << coos[coo][i] << " ";
		// 	std::cout << std::endl;
		// }
		// std::cout << std::endl;

		// for (int coo=0; coo<3; coo++){
		// 	std::vector<double> reordered_coo = reorder(coos[coo],node_order);
		// 	for (size_t i=0; i<50000; i++)
		// 		std::cout << reordered_coo[i] << " ";
		// 	std::cout << std::endl;
		// }
		// std::cout << std::endl;

		// return 0;

		///////////////////////////////////////////////////////////////


		if (config["max_block_id"]>=0 and block_id++>=config["max_block_id"]) break;
	}

	bpReader.EndStep();	// end logical step
	bpWriter.EndStep();	// end logical step
	}

	std::cout << "Total time = " << std::setprecision(3) << total_time << std::endl;

	bpReader.Close();	// close engine
	bpWriter.Close();	// close engine

#if ADIOS2_USE_MPI
	MPI_Finalize();
#endif
	return 0;
}
