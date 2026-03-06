#include <cmath>
#include <cstddef>

#include <algorithm>
#include <array>
#include <iostream>
#include <fstream>
#include <vector>
#include <numeric>

#include "mgard/compress.hpp"

#include "zfp.h"

#include <adios2.h>

#include <nlohmann/json.hpp>


using json = nlohmann::json;
using BYTE = unsigned char;


#define DTYPE double



unsigned char * zfp_compress(double * array, double tolerance, size_t r1, size_t r2, size_t r3, size_t *out_size){
    zfp_type type;     /* array scalar type */
    zfp_field* field;  /* array meta data */
    zfp_stream* zfp;   /* compressed stream */
    void* buffer;      /* storage for compressed stream */
    size_t bufsize;    /* byte size of compressed buffer */
    bitstream* stream; /* bit stream to write to or read from */
    size_t zfpsize;    /* byte size of compressed stream */

    /* allocate meta data for the 3D array a[nz][ny][nx] */
    type = zfp_type_double;
    if ((r1==1) && (r2==1)) {
        field = zfp_field_1d(array, type, r3);
    } else if (r1==1) {
        field = zfp_field_2d(array, type, r2, r3);
    }
    else {
        field = zfp_field_3d(array, type, r1, r2, r3);
    }
    /* allocate meta data for a compressed stream */
    zfp = zfp_stream_open(NULL);

    /* set compression mode and parameters via one of three functions */
    /*  zfp_stream_set_rate(zfp, rate, type, 3, 0); */
    /*  zfp_stream_set_precision(zfp, precision); */
    zfp_stream_set_accuracy(zfp, tolerance);

    /* allocate buffer for compressed data */
    bufsize = zfp_stream_maximum_size(zfp, field);
    buffer = malloc(bufsize);

    /* associate bit stream with allocated buffer */
    stream = stream_open(buffer, bufsize);
    zfp_stream_set_bit_stream(zfp, stream);
    zfp_stream_rewind(zfp);

    // std::cout << field->nx << " " << field->ny << " " << field->nz << " " << r3 << std::endl;

    zfpsize = zfp_compress(zfp, field);

    zfp_field_free(field);
    zfp_stream_close(zfp);
    stream_close(stream);
    *out_size = zfpsize;
    return (unsigned char *)buffer;
}


double * zfp_decompress(unsigned char * comp_data, double tolerance, size_t buffer_size, size_t r1, size_t r2, size_t r3){
    zfp_type type;     /* array scalar type */
    zfp_field* field;  /* array meta data */
    zfp_stream* zfp;   /* compressed stream */
    void* buffer;      /* storage for compressed stream */
    size_t bufsize;    /* byte size of compressed buffer */
    bitstream* stream; /* bit stream to write to or read from */

    /* allocate meta data for the 3D array a[nz][ny][nx] */
    double * array = (double *) malloc(r1 * r2 * r3 * sizeof(double));
    type = zfp_type_double;
    if ((r1==1) &&(r2==1)) {
        field = zfp_field_1d(array, type, r3);
    } else if (r1==1) {
        field = zfp_field_2d(array, type, r3, r2);
    }
    else
        field = zfp_field_3d(array, type, r3, r2, r1);

    /* allocate meta data for a compressed stream */
    zfp = zfp_stream_open(NULL);

    /* set compression mode and parameters via one of three functions */
    /*  zfp_stream_set_rate(zfp, rate, type, 3, 0); */
    /*  zfp_stream_set_precision(zfp, precision); */
    zfp_stream_set_accuracy(zfp, tolerance);

    /* allocate buffer for compressed data */
    bufsize = zfp_stream_maximum_size(zfp, field);
    // buffer = malloc(bufsize);
    buffer = (void *) comp_data;
    bufsize = buffer_size;

    /* associate bit stream with allocated buffer */
    stream = stream_open(buffer, bufsize);
    zfp_stream_set_bit_stream(zfp, stream);
    zfp_stream_rewind(zfp);

    zfp_decompress(zfp, field);

    zfp_field_free(field);
    zfp_stream_close(zfp);
    stream_close(stream);
    return array;
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

	// adios2::IO writer2_io = adios.DeclareIO("BPWriter2");
	// writer2_io.SetEngine("BPFile");
	// adios2::Engine bpWriter2 = writer2_io.Open(argv[4], adios2::Mode::Write);


	///////////////////////////////////////////////////////////////////////////
	// compression config

	const DTYPE s = config["s"];
	const DTYPE rel_tol = config["rel_tol"];


	///////////////////////////////////////////////////////////////////////////
	// define output ADIOS variables


	adios2::Variable<int64_t> out_adios_connectivity = writer_io.DefineVariable<int64_t>(connectivity_name, {}, {}, {adios2::UnknownDim});

	std::vector<adios2::Variable<double>> out_adios_coos(coo_names.size());
	for (int i=0; i<coo_names.size(); i++)
		out_adios_coos[i] = writer_io.DefineVariable<double>(coo_names[i], {}, {}, {adios2::UnknownDim});

	std::vector<adios2::Variable<BYTE>> out_adios_vars(var_names.size());
	for (int i=0; i<var_names.size(); i++)
		out_adios_vars[i] = writer_io.DefineVariable<BYTE>(var_names[i], {}, {}, {adios2::UnknownDim});


	///////////////////////////////////////////////////////////////////////////


	double total_time = 0;

	std::vector<double> original_size(var_names.size(), 0.0);
	std::vector<double> compressed_size(var_names.size(), 0.0);
	std::vector<double> rel_errors(var_names.size(), 0.0);


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
	// compress variables

	size_t tot_N = 0;


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
		// MGARD mesh

		// number of dimensions
		size_t ndims = coos.size();

		// number of nodes
		size_t N = coos[0].size();
		tot_N += N;

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
			// cumul_dst += 1;
			mgard_x.push_back(cumul_dst);
		}
		for (int i=0; i<N; i++){
			// mgard_x[i] /= mgard_x[N-1];
			mgard_x[i] = float(i)/(N-1);
		}

		// wrap the information about the mesh into an `mgard::TensorMeshHierarchy`
		const mgard::TensorMeshHierarchy<1, DTYPE> hierarchy({N}, mgard_coos);
		// const mgard::TensorMeshHierarchy<1, DTYPE> hierarchy({N});


		//////////////////////////////////////////////////////////////////////////////////
		// save compressed variables

		out_adios_connectivity.SetSelection(adios2::Box<adios2::Dims>({}, {ElementConnectivity.size()}));
		bpWriter.Put<int64_t>(out_adios_connectivity, ElementConnectivity.data(), adios2::Mode::Sync);

		for (int i=0; i<coo_names.size(); i++){
			out_adios_coos[i].SetSelection(adios2::Box<adios2::Dims>({}, {coos[i].size()}));
			bpWriter.Put<DTYPE>(out_adios_coos[i], coos[i].data(), adios2::Mode::Sync);
		}


		double avg_err = 0;

		for (int i=0; i<var_names.size(); i++){
			// find absolute tolerance
			DTYPE mag_v = 0;
			for (size_t k=0; k<N; k++)
				mag_v += vars[i][k] * vars[i][k] / N;
			// for (size_t k=1; k<N; k++)
			// 	mag_v += (vars[i][k-1] * vars[i][k-1] + vars[i][k] * vars[i][k]) * (mgard_x[k]-mgard_x[k-1]) / 2;
			double L2_norm = std::sqrt(mag_v);

			// double rel_tol_a = rel_tol;
			// double rel_tol_b = 20*rel_tol;
			// double rel_tol_i = 0.5*(rel_tol_a+rel_tol_b);
			// double rel_err;
			// for (int it=0; it<10; it++){
			// 	DTYPE abs_tol = rel_tol_i * L2_norm;

			// 	// // check for nan
			// 	// for (int j=1; j<N; j++)
			// 	// 	if(std::isnan(vars[i][j])) vars[i][j]=0;

			// 	// compress variable
			// 	const mgard::CompressedDataset<1, DTYPE> compressed = mgard::compress(hierarchy, vars[i].data(), s, abs_tol);

			// 	original_size[i]   += static_cast<DTYPE>(N*sizeof(DTYPE));
			// 	compressed_size[i] += compressed.size();

			// 	////////////////////////////////////////////////////////////
			// 	// Compression error
			// 	// if (var_names[i]=="/hpMusic_base/hpMusic_Zone/FlowSolution/U_aver"){

			// 	DTYPE *const u_copy = new DTYPE[N];
			// 	std::copy(vars[i].data(), vars[i].data() + N, u_copy);

			// 	const mgard::DecompressedDataset<1, DTYPE> decompressed = mgard::decompress(compressed);
			// 	DTYPE *const error = new DTYPE[N];
			// 	// DTYPE *const Linf  = new DTYPE[N];
			// 	// The `data` member function returns a pointer to the decompressed dataset.
			// 	std::transform(u_copy, u_copy + N, decompressed.data(), error, std::minus<DTYPE>());
			// 	delete[] u_copy;

			// 	DTYPE Linf;
			// 	double maxLinf = -10;
			// 	for (int j=0; j<N; j++){
			// 		Linf = std::fabs(error[j]);
			// 		if (Linf>maxLinf) maxLinf = Linf;
			// 	}

			// 	DTYPE *const shuffled = new DTYPE[N];
			// 	mgard::shuffle(hierarchy, error, shuffled);
			// 	delete[] error;
			// 	double L2_err = mgard::norm(hierarchy, shuffled, s);
			// 	rel_err = L2_err / L2_norm;

			// 	// std::cout << std::endl << var_names[i] << ": " << std::endl;
			// 	// std::cout << "rel. err. tolerance: " << rel_tol_i << std::endl
			// 	//           << "abs. err. tolerance: " << abs_tol << std::endl
			// 	//           << "achieved   L2 error: " << L2_err << std::endl
			// 	//           << "achieved  rel error: " << rel_err << std::endl
			// 	//           << "achieved Linf error: " << maxLinf
			// 	//           << std::endl;
			// 	// delete[] shuffled;

			// 	// // `compressed` contains the compressed data buffer. We can query its size in bytes with the `size` member function.
			// 	// std::cout << "  compression ratio: "
			// 	//           << static_cast<DTYPE>(N*sizeof(DTYPE)) / compressed.size()
			// 	//           << std::endl << std::endl;
			// 	////////////////////////////////////////////////////////////


			// 	if (rel_err<rel_tol){
			// 		if (std::abs(rel_err-rel_tol)<0.01*rel_tol) break;
			// 		rel_tol_a = rel_tol_i;
			// 		rel_tol_i = 0.5 * (rel_tol_a + rel_tol_b);
			// 	}else{
			// 		rel_tol_b = rel_tol_i;
			// 		rel_tol_i = 0.5 * (rel_tol_a + rel_tol_b);;
			// 	}
			// }

			double rel_err = (double)config["rel_tol"];
			DTYPE abs_tol = rel_err * L2_norm;

			// std::cout << var_names[i] << std::endl;
			// std::cout << L2_norm << std::endl;
			// std::cout << rel_err << std::endl;
			// std::cout << abs_tol << std::endl;
			// return 0;

			// compress variable
			// const mgard::CompressedDataset<1, DTYPE> compressed = mgard::compress(hierarchy, vars[i].data(), s, abs_tol);

			size_t out_size;
			unsigned char *compressed_data = zfp_compress(vars[i].data(), abs_tol, 1, 1, N, &out_size);
			double* decompressed_data = zfp_decompress(compressed_data, abs_tol, out_size, 1, 1, N);

			original_size[i]   += static_cast<DTYPE>(N*sizeof(DTYPE));
			compressed_size[i] += out_size;
			rel_errors[i] += rel_err / ((int)config["max_block_id"]+1); // / nblocks

			out_adios_vars[i].SetSelection(adios2::Box<adios2::Dims>({}, {out_size}));
			bpWriter.Put<BYTE>(out_adios_vars[i], (BYTE*)compressed_data, adios2::Mode::Sync);

			mag_v = 0;
			for (size_t k=0; k<N; k++){
				mag_v += (*((DTYPE*)decompressed_data+k)-vars[i][k]) * ( *((DTYPE*)decompressed_data+k)-vars[i][k]) / N;
			}
			double L2_norm_err = std::sqrt(mag_v);

			avg_err += L2_norm_err / (L2_norm+1.e-5) / vars.size();
		}
		std::cout << "\nerror " << avg_err << std::endl;
		bpWriter.PerformPuts();
		// bpWriter2.PerformPuts();


		///////////////////////////////////////////////////////////////

		auto end = std::chrono::high_resolution_clock::now();
		auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

		total_time += (double)duration.count() / 1e6;

		std::cout << std::setprecision(3) << ", time = " << (double)duration.count() / 1e6 << " sec" << std::endl;

		///////////////////////////////////////////////////////////////


		if (config["max_block_id"]>=0 and block_id++>=config["max_block_id"]) break;
	}


	// print compression ratios
	std::cout << " var_name:  "; for (int i=0; i<var_names.size(); i++) std::cout << var_names[i]     << " "; std::cout << std::endl;
	std::cout << " rel_error: "; for (int i=0; i<var_names.size(); i++) std::cout << rel_errors[i]    << " "; std::cout << std::endl;
	std::cout << " orig_size: "; for (int i=0; i<var_names.size(); i++) std::cout << original_size[i] << " ";
	std::cout << std::accumulate(original_size.begin(), original_size.end(), decltype(original_size)::value_type(0));
	std::cout << std::endl;
	std::cout << " comp_size: "; for (int i=0; i<var_names.size(); i++) std::cout << compressed_size[i] << " ";
	std::cout << std::accumulate(compressed_size.begin(), compressed_size.end(), decltype(compressed_size)::value_type(0));
	std::cout << std::endl;
	std::cout << " CR:        "; for (int i=0; i<var_names.size(); i++) std::cout << original_size[i] / compressed_size[i] << " ";
	std::cout << std::accumulate(original_size.begin(), original_size.end(), decltype(original_size)::value_type(0)) / std::accumulate(compressed_size.begin(), compressed_size.end(), decltype(compressed_size)::value_type(0)) << " ";
	std::cout << std::endl;


	std::cout << "\ntot_N " << tot_N << std::endl;

	// bpReader.EndStep(); // end logical step
	// bpReader.Close();  // close engine
	// bpWriter.Close();  // close engine
	// bpWriter2.Close();  // close engine

	bpReader.EndStep();	// end logical step
	bpWriter.EndStep();	// end logical step
	}


	std::cout << "\nTotal CR = " << std::setprecision(3)
			<< std::accumulate(original_size.begin(), original_size.end(), decltype(original_size)::value_type(0)) / std::accumulate(compressed_size.begin(), compressed_size.end(), decltype(compressed_size)::value_type(0))
			<< std::endl;
	std::cout << "Total time = " << std::setprecision(3) << total_time << std::endl;


	bpReader.Close();  // close engine
	bpWriter.Close();  // close engine

#if ADIOS2_USE_MPI
	MPI_Finalize();
#endif
	return 0;
}
