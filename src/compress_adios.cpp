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
	adios2::ADIOS adios;
	adios2::IO reader_io = adios.DeclareIO("BPReader");
	adios2::IO writer_io = adios.DeclareIO("BPWriter");
	reader_io.SetEngine("BPFile");
	writer_io.SetEngine("BPFile");
	adios2::Engine bpReader = reader_io.Open(argv[2], adios2::Mode::Read);
	adios2::Engine bpWriter = writer_io.Open(argv[3], adios2::Mode::Write);


	adios2::StepStatus read_status = bpReader.BeginStep(adios2::StepMode::Read, -1.0f);


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
	// define output ADIOS variables

	const DTYPE s = config["s"];
	const DTYPE rel_tol = config["rel_tol"];

	// ADIOS compression
	adios2::Operator op = adios.DefineOperator("mgardplus", "mgardplus");

	// adios2::Variable<int64_t> out_adios_connectivity = writer_io.DefineVariable<int64_t>(connectivity_name, {}, {}, {adios2::UnknownDim});

	// std::vector<adios2::Variable<double>> out_adios_coos(coo_names.size());
	// for (int i=0; i<coo_names.size(); i++)
	// 	out_adios_coos[i] = writer_io.DefineVariable<double>(coo_names[i], {}, {}, {adios2::UnknownDim});

	std::vector<adios2::Variable<DTYPE>> out_adios_vars(var_names.size());
	for (int i=0; i<var_names.size(); i++){
		// double abs_tol = rel_tol * (adios_vars[i].Max()-adios_vars[i].Min());
		// abs_tol = std::max(1.e-6, rel_tol * (adios_vars[i].Max()-adios_vars[i].Min()));
		out_adios_vars[i] = writer_io.DefineVariable<DTYPE>(var_names[i], {}, {}, {adios2::UnknownDim});
		// out_adios_vars[i].AddOperation(op, {{"accuracy", std::to_string(abs_tol)}, {"mode", "ABS"}});
		// std::cout << abs_tol << std::endl;
	}

	///////////////////////////////////////////////////////////////////////////
	// compress variables


	// number of blocks/subdomains
	auto blocks = bpReader.BlocksInfo(adios_connectivity, 0);
	size_t nblocks = blocks.size();
	int block_id = 0;
	for (auto &block : blocks){
		std::cout << "blockID = " << block.BlockID << " out of " << nblocks << std::endl;

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

		// out_adios_connectivity.SetSelection(adios2::Box<adios2::Dims>({}, {ElementConnectivity.size()}));
		// bpWriter.Put<int64_t>(out_adios_connectivity, ElementConnectivity.data(), adios2::Mode::Sync);

		// for (int i=0; i<coo_names.size(); i++){
		// 	out_adios_coos[i].SetSelection(adios2::Box<adios2::Dims>({}, {coos[i].size()}));
		// 	bpWriter.Put<DTYPE>(out_adios_coos[i], coos[i].data(), adios2::Mode::Sync);
		// }

		for (int i=0; i<var_names.size(); i++){
			// find absolute tolerance
			DTYPE mag_v = 0;
			for (size_t k=0; k<vars[i].size(); k++)
				mag_v += vars[i][k] * vars[i][k] / vars[i].size();
			DTYPE abs_tol = rel_tol * std::sqrt(mag_v);

			out_adios_vars[i].AddOperation(op, {{"accuracy", std::to_string(abs_tol)}, {"mode", "ABS"}});

			out_adios_vars[i].SetSelection(adios2::Box<adios2::Dims>({}, {vars[i].size()}));
			bpWriter.Put<DTYPE>(out_adios_vars[i], vars[i].data(), adios2::Mode::Sync);

			out_adios_vars[i].RemoveOperations();
		}
		bpWriter.PerformPuts();


	if (config["max_block_id"]>=0 and block_id++>=config["max_block_id"]) break;

	}

	bpReader.EndStep(); // end logical step
	bpReader.Close();  // close engine
	bpWriter.Close();  // close engine

}
