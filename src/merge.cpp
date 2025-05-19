/*
 * merge.cpp :
 *
 *  Created on: October 2, 2024
 *      Author: Viktor Reshniak
 *
 *	Merge DG data to make it continuous, this results in reindexing DOFs
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
		std::cerr << "Usage: merge config.json orig.bp merged.bp, merge got " << argc-1 << " parameters" << std::endl;
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


	size_t total_nodes = 0;
	size_t total_unique_nodes = 0;
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


	// number of blocks/subdomains
	auto blocks = bpReader.BlocksInfo(adios_connectivity, 0);
	size_t nblocks = blocks.size();

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
		// reindexing

		// Calculate inverse index map: unique linear index of the node -> all indices of this (duplicate) node
		auto inv_map = inverse_index_map(coos);

		// Convert inverse index map to the forward map: index of the node -> unique index of the node
		auto ind_map = index_map_from_inverse_map(inv_map, coos[0].size());

		// // Convert inverse map to serial vector for saving in .bp file
		// auto serial_inv_map = serialize_inv_map(inv_map, coos[0].size());

		// reindex dofs to unique indices
		for(size_t i = 0; i < ElementConnectivity.size(); i++){
			ElementConnectivity[i] = ind_map[ElementConnectivity[i]];
		}

		// // print Element Connectivity
		// {
		// std::cout << "Element Connectivity:" << std::endl;
		// int nnodes = config["n_nodes_in_element"];
		// for (int i=0; i<10*nnodes; i+=nnodes)
		// {
		// 	for (int j=0; j<nnodes; j++) std::cout << " " << ElementConnectivity[i+j];
		// 	std::cout << std::endl;
		// }
		// }
		// return 0;

		// // print variable at duplicate nodes
		// for (int vari=0; vari<vars.size(); vari++){
		// 	for (int i=0; i<inv_map[0].size(); i++){
		// 		std::cout << vars[vari][inv_map[0][i]] << " ";
		// 	}
		// 	std::cout << std::endl;
		// }
		// return 0;


		//////////////////////////////////////////////////////////////////////////////////
		// merge variables, DG to CG

		std::vector<std::vector<double>> merged_coos(coo_names.size());
		std::vector<std::vector<double>> merged_vars(var_names.size());

		for (int i=0; i<coo_names.size(); i++) merged_coos[i] = merge_values(inv_map, coos[i]);
		for (int i=0; i<var_names.size(); i++) merged_vars[i] = merge_values(inv_map, vars[i]);


		//////////////////////////////////////////////////////////////////////////////////
		// save merged variables

		// set block/subregion for the variables
		out_adios_connectivity.SetSelection(adios2::Box<adios2::Dims>({}, {ElementConnectivity.size()}));
		bpWriter.Put<int64_t>(out_adios_connectivity, ElementConnectivity.data(), adios2::Mode::Sync);
		for (int i=0; i<coo_names.size(); i++){
			out_adios_coos[i].SetSelection(adios2::Box<adios2::Dims>({}, {merged_coos[i].size()}));
			bpWriter.Put<double>(out_adios_coos[i], merged_coos[i].data(), adios2::Mode::Sync);
		}
		for (int i=0; i<var_names.size(); i++){
			out_adios_vars[i].SetSelection(adios2::Box<adios2::Dims>({}, {merged_vars[i].size()}));
			bpWriter.Put<double>(out_adios_vars[i], merged_vars[i].data(), adios2::Mode::Sync);
		}

		// inv_map_out.SetSelection(adios2::Box<adios2::Dims>({}, {serial_inv_map.size()})); bpWriter.Put<size_t>(inv_map_out, serial_inv_map.data(), adios2::Mode::Sync);
		bpWriter.PerformPuts();


		///////////////////////////////////////////////////////////////

		auto end = std::chrono::high_resolution_clock::now();
		auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

		total_nodes += coos[0].size();
		total_unique_nodes += inv_map.size();
		total_time += (double)duration.count() / 1e6;

		std::cout << "\tNum. orig./uniq. nodes = " << coos[0].size() << "/" << inv_map.size()
				<< std::setprecision(3)
				<< ", redundancy = " <<  (double)coos[0].size()/(double)inv_map.size()
				<< ", time = " << (double)duration.count() / 1e6 << " sec" << std::endl;

		// // print inverse map
		// {
		// 	std::cout << "\tInverse map:" << std::endl;
		// 	for (size_t i=0; i<10; i++)
		// 	{
		// 		std::cout << "\t\t" << i << ": [";
		// 		for (const auto & n : inv_map[i]) std::cout << " " << n;
		// 		std::cout << " ]" << std::endl;
		// 	}
		// }

		// // print forward map
		// {
		// 	std::cout << "\tForward map:" << std::endl;
		// 	for (size_t i=0; i<10; i++)
		// 		std::cout << "\t\t" << i << ": " << ind_map[i] << std::endl;
		// }

		///////////////////////////////////////////////////////////////

		// for (size_t i=0; i<10000; i++)
		// 	std::cout << vars[1][i] << " ";
		// std::cout << std::endl;

		// std::vector<double> unmerged_var(vars[1].size());
		// for (size_t i=0; i<merged_vars[1].size(); i++)
		// 	for (size_t j=0; j<inv_map[i].size(); j++)
		// 		unmerged_var[inv_map[i][j]] = merged_vars[1][i];
		// for (size_t i=0; i<10000; i++)
		// 	std::cout << unmerged_var[i] << " ";
		// std::cout << std::endl;

		///////////////////////////////////////////////////////////////

		if (config["max_block_id"]>=0 and block.BlockID>=config["max_block_id"]) break;
	}

	bpReader.EndStep();	// end logical step
	bpWriter.EndStep();	// end logical step
	}

	std::cout << "Total num. orig. nodes = " << total_nodes << std::endl;
	std::cout << "Total num. uniq. nodes = " << total_unique_nodes << std::endl;
	std::cout << "Total redundancy = " << (double)total_nodes/(double)total_unique_nodes << std::endl;
	std::cout << "Total time = " << std::setprecision(3) << total_time << std::endl;

	bpReader.Close();	// close engine
	bpWriter.Close();	// close engine


#if ADIOS2_USE_MPI
	MPI_Finalize();
#endif
	return 0;
}
