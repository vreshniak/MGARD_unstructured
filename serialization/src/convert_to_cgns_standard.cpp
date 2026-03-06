/*
 * merge.cpp :
 *
 *  Created on: November 12, 2025
 *      Author: Viktor Reshniak
 *
 *	Convert Connectivity to CGNS standard
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
		std::cerr << "Usage: merge config.json orig.bp merged.bp, merge got " << argc-1 << " parameters" << std::endl;
		return -1;
	}
	std::string config_file_name = argv[1];
	std::string input_file_name  = argv[2];
	std::string output_file_name = argv[3];

	std::ifstream f(config_file_name);
	json config = json::parse(f);


	// variable names
	std::vector<std::string> coo_names = config["coo_names"];
	std::vector<std::string> var_names = config["var_names"];
	std::string connectivity_name      = config["connectivity"];

	int nodes_per_el = config["n_nodes_in_element"];

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
	adios2::Variable<int64_t> orig_adios_connectivity = writer_io.DefineVariable<int64_t>(connectivity_name, {}, {}, {adios2::UnknownDim});
	adios2::Variable<int64_t> out_adios_connectivity = writer_io.DefineVariable<int64_t>("/ExplicitElem/Connectivity", {}, {}, {adios2::UnknownDim});
	adios2::Variable<int32_t> out_adios_elemnumnodes = writer_io.DefineVariable<int32_t>("/ExplicitElem/NumNodes", {}, {}, {adios2::UnknownDim});
	adios2::Variable<uint8_t> out_adios_elemtype     = writer_io.DefineVariable<uint8_t>("/ExplicitElem/Types", {}, {}, {adios2::UnknownDim});

	std::vector<adios2::Variable<double>> out_adios_coos(coo_names.size());
	for (int i=0; i<coo_names.size(); i++) out_adios_coos[i] = writer_io.DefineVariable<double>(coo_names[i], {}, {}, {adios2::UnknownDim});

	std::vector<adios2::Variable<double>> out_adios_vars(var_names.size());
	for (int i=0; i<var_names.size(); i++) out_adios_vars[i] = writer_io.DefineVariable<double>(var_names[i], {}, {}, {adios2::UnknownDim});

	///////////////////////////////////////////////////////////////////////////


	size_t total_nodes{};
	size_t total_unique_nodes{};
	double total_time{};


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
		for (int i=0; i<coo_names.size(); i++) adios_coos[i] = reader_io.InquireVariable<double>(coo_names[i]);

		std::vector<adios2::Variable<double>> adios_vars(var_names.size());
		for (int i=0; i<var_names.size(); i++) adios_vars[i] = reader_io.InquireVariable<double>(var_names[i]);

		///////////////////////////////////////////////////////////////////////////


		// number of blocks/subdomains
		auto blocks = bpReader.BlocksInfo(adios_connectivity, 0);
		size_t nblocks = blocks.size();
		if (config["max_block_id"]>0){
			nblocks = config["max_block_id"];
			nblocks++;
		}

		// Build per-rank block indices (round-robin split)
		std::vector<size_t> myBlockIdx;
		myBlockIdx.reserve((nblocks + nranks - 1) / nranks);
		for (size_t i = rank; i < nblocks; i += nranks)
			myBlockIdx.push_back(i);


		for (size_t bi : myBlockIdx){
			const auto &block = blocks[bi];

			auto start = std::chrono::high_resolution_clock::now();
			if (verbose && rank == 0) std::cout << "step = " << step << ", block " << std::setw(2) << block.BlockID+1 << " out of " << nblocks;

			//////////////////////////////////////////////////////////////////////////////////
			// read from .bp file

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

			size_t num_elements = ElementConnectivity.size() / nodes_per_el;

			std::vector<int64_t> out_ElementConnectivity(ElementConnectivity.size());
			std::vector<int32_t> out_ElementNumNodes(num_elements);
			std::vector<uint8_t> out_ElementType(num_elements);

			int64_t out_ElementConnectivity_size{0};

			// first element
			out_ElementConnectivity[0] = ElementConnectivity[0];
			auto prev_id{ElementConnectivity[0]};
			int uniq_nodes_per_el = 1;
			for(int i = 1; i < nodes_per_el; i++){
				auto curr_id{ElementConnectivity[i]};
				if (curr_id != prev_id){
					out_ElementConnectivity[uniq_nodes_per_el] = curr_id;
					uniq_nodes_per_el++;
				}
				prev_id = curr_id;
			}
			out_ElementNumNodes[0] = uniq_nodes_per_el;
			switch (uniq_nodes_per_el) {
				case 4:
					out_ElementType[0] = 10;
					break;
				case 8:
					out_ElementType[0] = 12;
					break;
				case 6:
					out_ElementType[0] = 13;
					break;
				case 5:
					out_ElementType[0] = 14;
					break;
			}
			out_ElementConnectivity_size += uniq_nodes_per_el;


			// go through remaining elements
			for(auto k = nodes_per_el, elem_id = 1; k < ElementConnectivity.size(); k+=nodes_per_el, elem_id++){
				out_ElementConnectivity[out_ElementConnectivity_size] = ElementConnectivity[k];
				auto prev_id{ElementConnectivity[k]};
				uniq_nodes_per_el = 1;
				for(int i = 1; i < nodes_per_el; i++){
					auto curr_id{ElementConnectivity[k+i]};
					if (curr_id != prev_id){
						out_ElementConnectivity[out_ElementConnectivity_size+uniq_nodes_per_el] = curr_id;
						uniq_nodes_per_el++;
					}
					prev_id = curr_id;
				}
				out_ElementNumNodes[elem_id] = uniq_nodes_per_el;
				switch (uniq_nodes_per_el) {
					case 4:
						out_ElementType[elem_id] = 10;
						break;
					case 8:
						out_ElementType[elem_id] = 12;
						break;
					case 6:
						out_ElementType[elem_id] = 13;
						break;
					case 5:
						out_ElementType[elem_id] = 14;
						break;
				}
				out_ElementConnectivity_size += uniq_nodes_per_el;
			}


			//////////////////////////////////////////////////////////////////////////////////
			// save merged variables

			orig_adios_connectivity.SetSelection(adios2::Box<adios2::Dims>({}, {ElementConnectivity.size()}));
			bpWriter.Put<int64_t>(orig_adios_connectivity, ElementConnectivity.data(), adios2::Mode::Sync);

			out_adios_connectivity.SetSelection(adios2::Box<adios2::Dims>({}, {out_ElementConnectivity_size}));
			bpWriter.Put<int64_t>(out_adios_connectivity, out_ElementConnectivity.data(), adios2::Mode::Sync);

			out_adios_elemnumnodes.SetSelection(adios2::Box<adios2::Dims>({}, {out_ElementNumNodes.size()}));
			bpWriter.Put<int32_t>(out_adios_elemnumnodes, out_ElementNumNodes.data(), adios2::Mode::Sync);

			out_adios_elemtype.SetSelection(adios2::Box<adios2::Dims>({}, {out_ElementType.size()}));
			bpWriter.Put<uint8_t>(out_adios_elemtype, out_ElementType.data(), adios2::Mode::Sync);

			for (int i=0; i<coo_names.size(); i++){
				out_adios_coos[i].SetSelection(adios2::Box<adios2::Dims>({}, {coos[i].size()}));
				bpWriter.Put<double>(out_adios_coos[i], coos[i].data(), adios2::Mode::Sync);
			}
			for (int i=0; i<var_names.size(); i++){
				out_adios_vars[i].SetSelection(adios2::Box<adios2::Dims>({}, {vars[i].size()}));
				bpWriter.Put<double>(out_adios_vars[i], vars[i].data(), adios2::Mode::Sync);
			}

			bpWriter.PerformPuts();


			///////////////////////////////////////////////////////////////

			auto end = std::chrono::high_resolution_clock::now();
			auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

			// total_nodes += coos[0].size();
			// total_unique_nodes += inv_map.size();
			total_time += (double)duration.count() / 1e6;

			if (verbose && rank == 0)
			std::cout  << ", num_elems = " << num_elements << ", time = " << std::setprecision(2) << (double)duration.count() / 1e6 << " sec" << std::endl;

			///////////////////////////////////////////////////////////////

			if (config["max_block_id"]>=0 and block.BlockID>=config["max_block_id"]) break;
		} // end for myBlockIdx


		///////////////////////////////////////////////////////////////

		bpReader.EndStep();	// end logical step
		bpWriter.EndStep();	// end logical step

	} // end for time step

	if (verbose){
		std::cout << "Rank " << rank << std::endl;
		std::cout << "Total time = " << std::setprecision(3) << total_time << std::endl;
	}

	bpReader.Close();	// close engine
	bpWriter.Close();	// close engine


#if ADIOS2_USE_MPI
	MPI_Finalize();
#endif
	return 0;
}
