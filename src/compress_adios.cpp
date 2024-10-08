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


	///////////////////////////////////////////////////////////////////////////
	// compression config


	const DTYPE s = config["s"];
	const DTYPE rel_tol = config["rel_tol"];

	// ADIOS compression
	adios2::Operator op = adios.DefineOperator("mgardplus", "mgardplus");


	///////////////////////////////////////////////////////////////////////////
	// define output ADIOS variables


	adios2::Variable<int64_t> out_adios_connectivity = writer_io.DefineVariable<int64_t>(connectivity_name, {}, {}, {adios2::UnknownDim});

	std::vector<adios2::Variable<double>> out_adios_coos(coo_names.size());
	for (int i=0; i<coo_names.size(); i++)
		out_adios_coos[i] = writer_io.DefineVariable<double>(coo_names[i], {}, {}, {adios2::UnknownDim});

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

	adios2::Variable<int64_t> adios_connectivity = reader_io.InquireVariable<int64_t>(connectivity_name);

	std::vector<adios2::Variable<double>> adios_coos(coo_names.size());
	for (int i=0; i<coo_names.size(); i++)
		adios_coos[i] = reader_io.InquireVariable<double>(coo_names[i]);

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
	auto blocks = bpReader.BlocksInfo(adios_connectivity, 0);
	size_t nblocks = blocks.size();
	int block_id = 0;
	for (auto &block : blocks){
		auto start = std::chrono::high_resolution_clock::now();
		std::cout << "step = " << step << ", blockID = " << block.BlockID << "/" << nblocks;

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
		// save compressed variables

		out_adios_connectivity.SetSelection(adios2::Box<adios2::Dims>({}, {ElementConnectivity.size()}));
		bpWriter.Put<int64_t>(out_adios_connectivity, ElementConnectivity.data(), adios2::Mode::Sync);

		for (int i=0; i<coo_names.size(); i++){
			out_adios_coos[i].SetSelection(adios2::Box<adios2::Dims>({}, {coos[i].size()}));
			bpWriter.Put<DTYPE>(out_adios_coos[i], coos[i].data(), adios2::Mode::Sync);
		}

		for (int i=0; i<var_names.size(); i++){
			// find absolute tolerance
			DTYPE mag_v = 0;
			for (size_t k=0; k<vars[i].size(); k++)
				mag_v += vars[i][k] * vars[i][k] / vars[i].size();
			DTYPE L2_norm = std::sqrt(mag_v);
			DTYPE abs_tol = rel_tol * L2_norm;

			out_adios_vars[i].AddOperation(op, {{"accuracy", std::to_string(abs_tol)}, {"mode", "ABS"}});

			out_adios_vars[i].SetSelection(adios2::Box<adios2::Dims>({}, {vars[i].size()}));
			bpWriter.Put<DTYPE>(out_adios_vars[i], vars[i].data(), adios2::Mode::Sync);

			out_adios_vars[i].RemoveOperations();
		}
		bpWriter.PerformPuts();


		///////////////////////////////////////////////////////////////

		auto end = std::chrono::high_resolution_clock::now();
		auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

		total_time += (double)duration.count() / 1e6;

		std::cout << std::setprecision(3) << ", time = " << (double)duration.count() / 1e6 << " sec" << std::endl;


		///////////////////////////////////////////////////////////////

		if (config["max_block_id"]>=0 and block_id++>=config["max_block_id"]) break;
	}

	bpReader.EndStep();	// end logical step
	bpWriter.EndStep();	// end logical step
	}


	std::cout << "\nTotal time = " << std::setprecision(3) << total_time << std::endl;


	bpReader.Close();  // close engine
	bpWriter.Close();  // close engine

#if ADIOS2_USE_MPI
	MPI_Finalize();
#endif
	return 0;
}
