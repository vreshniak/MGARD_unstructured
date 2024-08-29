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
	for (int i=0; i<var_names.size(); i++)
		adios_vars[i] = reader_io.InquireVariable<double>(var_names[i]);


	///////////////////////////////////////////////////////////////////////////
	// define output ADIOS variables

	// adios2::Variable<int64_t> out_adios_connectivity = writer_io.DefineVariable<int64_t>(connectivity_name, {}, {}, {adios2::UnknownDim});

	// std::vector<adios2::Variable<double>> out_adios_coos(coo_names.size());
	// for (int i=0; i<coo_names.size(); i++)
	// 	out_adios_coos[i] = writer_io.DefineVariable<double>(coo_names[i], {}, {}, {adios2::UnknownDim});

	std::vector<adios2::Variable<BYTE>> out_adios_vars(var_names.size());
	for (int i=0; i<var_names.size(); i++)
		out_adios_vars[i] = writer_io.DefineVariable<BYTE>(var_names[i], {}, {}, {adios2::UnknownDim});


	///////////////////////////////////////////////////////////////////////////
	// compress variables


	const DTYPE s = config["s"];
	const DTYPE rel_tol = config["rel_tol"];


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
		// MGARD mesh

		// number of dimensions
		size_t ndims = coos.size();

		// number of nodes
		size_t N = coos[0].size();

		// Coordinate array
		std::array<std::vector<DTYPE>, 1> mgard_coos;
		std::vector<DTYPE> &mgard_x = mgard_coos.at(0);
		mgard_x.reserve(N);
		mgard_x.push_back(0.0);
		double cumul_dst = 0;
		for (int i=1; i<N; i++){
			double dst = 0;
			for (auto &coo : coos)
				dst += (coo[i]-coo[i-1]) * (coo[i]-coo[i-1]);
			cumul_dst += std::sqrt(dst) + 1.e-10;
			mgard_x.push_back(cumul_dst);
		}

		// wrap the information about the mesh into an `mgard::TensorMeshHierarchy`
		const mgard::TensorMeshHierarchy<1, DTYPE> hierarchy({N}, mgard_coos);


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
			for (size_t k=0; k<N; k++)
				mag_v += vars[i][k] * vars[i][k] / N;
			DTYPE abs_tol = rel_tol * std::sqrt(mag_v);

			// check for nan
			for (int j=1; j<N; j++)
				if(std::isnan(vars[i][j])) vars[i][j]=0;

			// compress variable
			const mgard::CompressedDataset<1, DTYPE> compressed = mgard::compress(hierarchy, vars[i].data(), s, abs_tol);

			////////////////////////////////////////////////////////////
			// Compression error
			if (var_names[i]=="FlowSolution/RHE"){
				DTYPE *const u_copy = new DTYPE[N];
				std::copy(vars[i].data(), vars[i].data() + N, u_copy);

				const mgard::DecompressedDataset<1, DTYPE> decompressed = mgard::decompress(compressed);
				DTYPE *const error = new DTYPE[N];
				// DTYPE *const Linf  = new DTYPE[N];
				// The `data` member function returns a pointer to the decompressed dataset.
				std::transform(u_copy, u_copy + N, decompressed.data(), error, std::minus<DTYPE>());
				delete[] u_copy;

				DTYPE Linf;
				double maxLinf = -10;
				for (int j=0; j<N; j++){
					Linf = std::fabs(error[j]);
					if (Linf>maxLinf) maxLinf = Linf;
				}

				DTYPE *const shuffled = new DTYPE[N];
				mgard::shuffle(hierarchy, error, shuffled);
				delete[] error;
				std::cout << std::endl << var_names[i] << ": " << std::endl;
				std::cout << "rel. err. tolerance: " << rel_tol << std::endl
				          << "abs. err. tolerance: " << abs_tol << std::endl
				          << "achieved   L2 error: " << mgard::norm(hierarchy, shuffled, s) << std::endl
				          << "achieved Linf error: " << maxLinf
				          << std::endl;
				delete[] shuffled;


				// `compressed` contains the compressed data buffer. We can query its size in bytes with the `size` member function.
				std::cout << "  compression ratio: "
				          << static_cast<DTYPE>(N*sizeof(DTYPE)) / compressed.size()
				          << std::endl << std::endl;

			}
			////////////////////////////////////////////////////////////

			out_adios_vars[i].SetSelection(adios2::Box<adios2::Dims>({}, {compressed.size()}));
			bpWriter.Put<BYTE>(out_adios_vars[i], (BYTE*)compressed.data(), adios2::Mode::Sync);

		}
		bpWriter.PerformPuts();


	if (config["max_block_id"]>=0 and block_id++>=config["max_block_id"]) break;

	}

	bpReader.EndStep(); // end logical step
	bpReader.Close();  // close engine
	bpWriter.Close();  // close engine

}
