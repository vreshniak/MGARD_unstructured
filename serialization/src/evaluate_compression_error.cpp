/*
 * evaluate_compression_error.cpp :
 *
 *  Created on: October 2, 2024
 *      Author: Viktor Reshniak
 *
 *	Reorder serialized DOFs
 */


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


bool verbose = true;


template <typename T>
double Linf(const std::vector<T> arr)
{
	return std::abs(*std::max_element(arr.begin(), arr.end(), [](double a, double b) {return std::abs(a) < std::abs(b);}));
}


template <typename T>
double L2(const std::vector<T> arr)
{
	double mag_v = 0;
	for (size_t k=0; k<arr.size(); k++)
		mag_v += std::pow(arr[k],2) / arr.size();
	return std::sqrt(mag_v);
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

	// read input
	if (argc<4){
		std::cerr << "Usage: evaluate_compression_error config.json orig.bp compressed.bp, `evaluate_compression_error` got " << argc-1 << " parameters" << std::endl;
		return -1;
	}
	std::string config_file_name = argv[1];
	std::string input1_file_name = argv[2];
	std::string input2_file_name = argv[3];

	// read config file
	std::ifstream f(config_file_name);
	json config = json::parse(f);

	// variable names
	std::vector<std::string> var_names = config["var_names"];

	int lenname = 0;
	for (auto &varname: var_names)
		if (varname.length()>lenname)
			lenname = varname.size();


	///////////////////////////////////////////////////////////////////////////
	// ADIOS readers/writers

	// setup ADIOS readers/writers
#if ADIOS2_USE_MPI
	adios2::ADIOS adios(MPI_COMM_WORLD);
#else
	adios2::ADIOS adios;
#endif
	adios2::IO reader1_io = adios.DeclareIO("BPReader1");
	adios2::IO reader2_io = adios.DeclareIO("BPReader2");
	reader1_io.SetEngine("BPFile");
	reader2_io.SetEngine("BPFile");
	adios2::Engine bpReader1 = reader1_io.Open(input1_file_name, adios2::Mode::Read);
	adios2::Engine bpReader2 = reader2_io.Open(input2_file_name, adios2::Mode::Read);


	///////////////////////////////////////////////////////////////////////////
	// define output ADIOS variables

	// std::vector<adios2::Variable<DTYPE>> out_adios_vars(var_names.size());
	// for (int i=0; i<var_names.size(); i++){
	// 	out_adios_vars[i] = writer_io.DefineVariable<DTYPE>(var_names[i], {}, {}, {adios2::UnknownDim});
	// }

	///////////////////////////////////////////////////////////////////////////
	// read - decompress - write

	// double total_time = 0;

	// time stepping
	while (true){
		adios2::StepStatus read_status1 = bpReader1.BeginStep(adios2::StepMode::Read, -1.0f);
		adios2::StepStatus read_status2 = bpReader2.BeginStep(adios2::StepMode::Read, -1.0f);
		if (read_status1!=adios2::StepStatus::OK) break;
		if (read_status2!=adios2::StepStatus::OK) break;

		size_t step = bpReader1.CurrentStep();

		///////////////////////////////////////////////////////////////////////////
		// define input ADIOS variables

		std::vector<adios2::Variable<double>> adios_vars1(var_names.size());
		for (int i=0; i<var_names.size(); i++)
			adios_vars1[i] = reader1_io.InquireVariable<double>(var_names[i]);

		std::vector<adios2::Variable<double>> adios_vars2(var_names.size());
		for (int i=0; i<var_names.size(); i++)
			adios_vars2[i] = reader2_io.InquireVariable<double>(var_names[i]);

		///////////////////////////////////////////////////////////////////////////
		// evaluate error for variables

		// number of blocks/subdomains
		auto blocks = bpReader1.BlocksInfo(adios_vars1[0], step); // blocks_info(name, step)
		size_t nblocks = blocks.size();
		for (auto &block : blocks){
			int block_id = block.BlockID;
			auto start = std::chrono::high_resolution_clock::now();
			if (verbose) std::cout << "step = " << step << ", blockID = " << block.BlockID << "/" << nblocks << std::endl;

			// allocate memory for reading variables
			std::vector<std::vector<double>> vars1(var_names.size());
			std::vector<std::vector<double>> vars2(var_names.size());

			for (auto &var : adios_vars1) var.SetBlockSelection(block.BlockID);
			for (auto &var : adios_vars2) var.SetBlockSelection(block.BlockID);

			// read variables
			for (int i=0; i<var_names.size(); i++){
				bpReader1.Get(adios_vars1[i], vars1[i], adios2::Mode::Sync);
				bpReader2.Get(adios_vars2[i], vars2[i], adios2::Mode::Sync);
			}
			bpReader1.PerformGets();
			bpReader2.PerformGets();

			for (int i=0; i<var_names.size(); i++){
				std::vector<double> err(vars1[i].size());
				std::transform(vars1[i].begin(), vars1[i].end(), vars2[i].begin(), err.begin(), [](double a, double b) { return a - b; });
				std::cout << std::left << std::setw(lenname) << var_names[i] << ": " << std::scientific << std::setprecision(2)
					<< "Linf= " << Linf<double>(err) / Linf<double>(vars1[i])
					<< "  L2= " << L2<double>(err) / L2<double>(vars1[i])
					<< std::endl;
					// << *std::max_element(err.begin(),      err.end(),      [](double a, double b) {return std::abs(a) < std::abs(b);}) << " / "
					// << *std::max_element(vars1[i].begin(), vars1[i].end(), [](double a, double b) {return std::abs(a) < std::abs(b);}) << std::endl;
			}


			///////////////////////////////////////////////////////////////

			if (config["max_block_id"]>=0 and block_id>=config["max_block_id"]) break;
		}

		bpReader1.EndStep();	// end logical step
		bpReader2.EndStep();	// end logical step

		if (step>=config["max_step"]) break;
	}


	bpReader1.Close();  // close engine
	bpReader2.Close();  // close engine

#if ADIOS2_USE_MPI
	MPI_Finalize();
#endif
	return 0;
}
