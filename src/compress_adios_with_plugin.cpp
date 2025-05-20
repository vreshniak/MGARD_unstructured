#include <cmath>
#include <cstddef>

#include <algorithm>
#include <array>
#include <iostream>
#include <fstream>
#include <vector>

#include "mgard/compress.hpp"

#include <adios2.h>

#include <nlohmann/json.hpp>


using json = nlohmann::json;
using BYTE = unsigned char;


#define DTYPE double


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
		std::cerr << "Usage: compress_adios_plugin config.json orig.bp compressed.bp, compress_adios_plugin got " << argc-1 << " parameters" << std::endl;
		return -1;
	}
	std::string config_file_name = argv[1];
	std::string input_file_name  = argv[2];
	std::string output_file_name = argv[3];

	std::ifstream f(config_file_name);
	json config = json::parse(f);


	// variable names
	std::vector<std::string> var_names = config["var_names"];


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
	// compression config


	const DTYPE s = config["s"];
	const DTYPE rel_tol = config["rel_tol"];

	// ADIOS compression
	// adios2::Operator op = adios.DefineOperator("mgardplus", "mgardplus");


	///////////////////////////////////////////////////////////////////////////
	// define output ADIOS variables

	std::vector<adios2::Variable<DTYPE>> out_adios_vars(var_names.size());
	for (int i=0; i<var_names.size(); i++){
		out_adios_vars[i] = writer_io.DefineVariable<DTYPE>(var_names[i], {}, {}, {adios2::UnknownDim});
	}

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

		std::vector<adios2::Variable<double>> adios_vars(var_names.size());
		std::vector<double> vars_min(var_names.size());
		std::vector<double> vars_max(var_names.size());
		for (int i=0; i<var_names.size(); i++){
			adios_vars[i] = reader_io.InquireVariable<double>(var_names[i]);
			vars_min[i] = adios_vars[i].Min();
			vars_max[i] = adios_vars[i].Max();
		}


		///////////////////////////////////////////////////////////////////////////
		// compress variables


		// number of blocks/subdomains
		auto blocks = bpReader.BlocksInfo(adios_vars[0], 0);
		size_t nblocks = blocks.size();
		int block_id = 0;
		for (auto &block : blocks){
			auto start = std::chrono::high_resolution_clock::now();
			std::cout << "step = " << step << ", blockID = " << block.BlockID << "/" << nblocks << std::endl;

			// allocate memory for reading variables
			std::vector<std::vector<double>> vars(var_names.size());

			// set block/subregion for the variables
			// adios_connectivity.SetBlockSelection(block.BlockID);
			for (auto &var : adios_vars) var.SetBlockSelection(block.BlockID);

			// read variables
			for (int i=0; i<var_names.size(); i++) bpReader.Get(adios_vars[i], vars[i], adios2::Mode::Sync);
			bpReader.PerformGets();


			//////////////////////////////////////////////////////////////////////////////////
			// save compressed variables

			for (int i=0; i<var_names.size(); i++){
				// find absolute tolerance
				DTYPE mag_v = 0;
				for (size_t k=0; k<vars[i].size(); k++)
					mag_v += vars[i][k] * vars[i][k] / vars[i].size();
				DTYPE L2_norm = std::sqrt(mag_v);
				DTYPE abs_tol = rel_tol * L2_norm;

				adios2::Params params;
				params["PluginName"] = "serialize_operator";
				params["PluginLibrary"] = "CompressMGARDSerializeOperator";
				params["meshfile"] = "/lustre/orion/cfd164/proj-shared/reshniakv/order.bp"; //argv[2];
				params["blockid"] = std::to_string(block_id);
				params["accuracy"] = std::to_string(abs_tol);
				out_adios_vars[i].AddOperation("plugin", params);
				// // out_adios_vars[i].AddOperation(op, {{"accuracy", std::to_string(abs_tol)}, {"mode", "ABS"}});

				out_adios_vars[i].SetSelection(adios2::Box<adios2::Dims>({}, {vars[i].size()}));
				bpWriter.Put<DTYPE>(out_adios_vars[i], vars[i].data(), adios2::Mode::Sync);

				out_adios_vars[i].RemoveOperations();
			}
			bpWriter.PerformPuts();


			///////////////////////////////////////////////////////////////

			// auto end = std::chrono::high_resolution_clock::now();
			// auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

			// total_time += (double)duration.count() / 1e6;

			// std::cout << std::setprecision(3) << ", time = " << (double)duration.count() / 1e6 << " sec" << std::endl;


			///////////////////////////////////////////////////////////////

			if (config["max_block_id"]>=0 and block_id++>=config["max_block_id"]) break;
		}

		bpReader.EndStep();	// end logical step
		bpWriter.EndStep();	// end logical step
	}


	// std::cout << "\nTotal time = " << std::setprecision(3) << total_time << std::endl;


	bpReader.Close();  // close engine
	bpWriter.Close();  // close engine

#if ADIOS2_USE_MPI
	MPI_Finalize();
#endif
	return 0;
}
