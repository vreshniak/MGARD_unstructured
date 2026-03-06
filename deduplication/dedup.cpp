/*
 * dedup.cpp :
 *
 *  Created on: March 5, 2026
 *      Author: Viktor Reshniak
 *
 *	Deduplicate DG data to make it continuous, this results in reindexing DOFs
 */


#include <adios2.h>
#include <vector>
#include <algorithm>    // std::find
#include <unordered_map>

#include <iostream>
#include <fstream>

#include <mpi.h>
#include <chrono>

#include <nlohmann/json.hpp>
using json = nlohmann::json;


bool verbose = true;
bool debug   = false;



/**
 * Vector hash functor for unordered_map with vector keys
 */
template <typename T>
struct VectorHash
{
	size_t operator()(const std::vector<T>& v) const noexcept
	{
		size_t h = 0;
		for (const auto& x : v) {
			size_t hx = std::hash<T>{}(x);
			// hash combine formula: h = h * 31 + hx; (or any other prime number)
			h ^= hx + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
		}
		return h;
	}
};


struct IndexMapResult
{
	std::vector<size_t> forward_map;
	size_t n_unique;
};


/**
 * Build a forward index map from coordinates: original index -> unique index
 *
 * Inputs:
 * ------
 *  coos - vector of x-y-z coordinate vectors
*/
template <typename Container>
IndexMapResult index_map_from_coos(const Container& coos)
{
	using scalar_t = typename Container::value_type::value_type;
	using node_t   = std::vector<scalar_t>;

	const size_t ndim = coos.size();
	const size_t n    = coos.empty() ? 0 : coos[0].size();

	// std::map<node_t, size_t> unique_index_of;
	std::unordered_map<node_t, size_t, VectorHash<scalar_t>> unique_index_of;
    unique_index_of.reserve(n);
	
	std::vector<size_t> forward_map(n);
	node_t curr(ndim);

	// go through all nodes and assign unique indices to them
	for (size_t i = 0; i < n; ++i){
		// ndim-dimensional node
		for (size_t d = 0; d < ndim; ++d) curr[d] = coos[d][i];

		auto [it, inserted] = unique_index_of.emplace(curr, unique_index_of.size());
		forward_map[i] = it->second;
	}

	size_t n_unique = unique_index_of.size();
	return {std::move(forward_map), n_unique};
}


/**
 * Merge values at duplicate nodes by averaging them
 *
 * Inputs:
 * ------
 *  forward_map - index of the node -> unique index of the node
 *  values - variable values at the original nodes (with duplicates)
 *  n_unique - number of unique nodes
 */
template <typename T>
std::vector<T> merge_values(const std::vector<size_t>& forward_map, const std::vector<T>& values, size_t n_unique)
{
	std::vector<T> sum(n_unique, T{});
	std::vector<size_t> count(n_unique, 0);

	const size_t n = values.size();
	// go through all nodes and sum values at duplicate nodes, also count the number of duplicates for averaging
	for (size_t i = 0; i < n; ++i){
		// unique index of the node
		const size_t u = forward_map[i];
		sum[u] += values[i];
		++count[u];
	}

	// average values at duplicate nodes
	for (size_t u = 0; u < n_unique; ++u){
		if (count[u] > 0) sum[u] /= static_cast<T>(count[u]);
	}
	return sum; // merged values
}



// lambda to get CGNS element type from number of unique nodes in the element
auto elem_type_from_num_nodes = [](int uniq_nodes) -> uint8_t {
	switch (uniq_nodes) {
		case 4: return 10;
		case 8: return 12;
		case 6: return 13;
		case 5: return 14;
		default: return 0;
	}
};



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
		std::cerr << "Usage: dedup config.json orig.bp merged.bp, dedup got " << argc-1 << " parameters" << std::endl;
		return -1;
	}
	std::string config_file_name = argv[1];
	std::string input_file_name  = argv[2];
	std::string output_file_name = argv[3];

	// parse config file
	std::ifstream f(config_file_name);
	json config = json::parse(f);

	// variable names
	std::string connectivity_name      = config["connectivity"];
	std::vector<std::string> coo_names = config["coo_names"];
	std::vector<std::string> var_names = config["var_names"];

	int nodes_per_el = config["n_nodes_in_element"];
	

	int64_t max_steps = config["max_step"];
	int64_t max_block_id = config["max_block_id"];
	size_t n_nodes_in_element = config["n_nodes_in_element"];

	bool convert_to_cgns_standard = config["convert_to_cgns_standard"];
	verbose = config["verbose"];

	///////////////////////////////////////////////////////////////////////////
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

	// adios2::Variable<int64_t> out_adios_connectivity = writer_io.DefineVariable<int64_t>(connectivity_name, {}, {}, {adios2::UnknownDim});
	adios2::Variable<int64_t> orig_adios_connectivity = writer_io.DefineVariable<int64_t>(connectivity_name, {}, {}, {adios2::UnknownDim});
	adios2::Variable<int64_t> out_adios_connectivity  = writer_io.DefineVariable<int64_t>("/ExplicitElem/Connectivity", {}, {}, {adios2::UnknownDim});
	adios2::Variable<int32_t> out_adios_elemnumnodes  = writer_io.DefineVariable<int32_t>("/ExplicitElem/NumNodes", {}, {}, {adios2::UnknownDim});
	adios2::Variable<uint8_t> out_adios_elemtype      = writer_io.DefineVariable<uint8_t>("/ExplicitElem/Types", {}, {}, {adios2::UnknownDim});

	std::vector<adios2::Variable<double>> out_adios_coos(coo_names.size());
	std::vector<adios2::Variable<double>> out_adios_vars(var_names.size());
	for (int i=0; i<coo_names.size(); i++) out_adios_coos[i] = writer_io.DefineVariable<double>(coo_names[i], {}, {}, {adios2::UnknownDim});
	for (int i=0; i<var_names.size(); i++) out_adios_vars[i] = writer_io.DefineVariable<double>(var_names[i], {}, {}, {adios2::UnknownDim});


	///////////////////////////////////////////////////////////////////////////
	// time stepping loop and block loop

	size_t total_nodes{};
	size_t total_unique_nodes{};
	double total_time{};

	// time stepping
	while (true){
		// begin step: this advances to the next time step and makes the variables available for reading
		adios2::StepStatus read_status = bpReader.BeginStep(adios2::StepMode::Read, -1.0f);
		size_t step = bpReader.CurrentStep();
		// break if we have read all steps, or if we have reached max_steps (when max_steps >= 0)
		if (read_status!=adios2::StepStatus::OK || (max_steps >= 0 && step > static_cast<size_t>(max_steps)))
			break;
		bpWriter.BeginStep();


		///////////////////////////////////////////////////////////////////////////
		// define input ADIOS variables
		// this needs to be done after BeginStep() because the variables are not available before that

		adios2::Variable<int64_t> adios_connectivity = reader_io.InquireVariable<int64_t>(connectivity_name);
		if (!adios_connectivity) throw std::runtime_error("Missing variable: " + connectivity_name);

		std::vector<adios2::Variable<double>> adios_coos(coo_names.size());
		std::vector<adios2::Variable<double>> adios_vars(var_names.size());
		for (int i=0; i<coo_names.size(); i++){
			adios_coos[i] = reader_io.InquireVariable<double>(coo_names[i]);
			if (!adios_coos[i]) throw std::runtime_error("Missing variable: " + coo_names[i]);
		}
		for (int i=0; i<var_names.size(); i++){
			adios_vars[i] = reader_io.InquireVariable<double>(var_names[i]);
			if (!adios_vars[i]) throw std::runtime_error("Missing variable: " + var_names[i]);
		}

		
		///////////////////////////////////////////////////////////////////////////
		// number of blocks/subdomains
		auto blocks = bpReader.BlocksInfo(adios_connectivity, 0);
		size_t nblocks = blocks.size();
		if (max_block_id >= 0){
			size_t capped_nblocks = static_cast<size_t>(max_block_id) + 1;
			if (capped_nblocks < nblocks)
				nblocks = capped_nblocks;
		}
		

		// loop over my blocks and read, deduplicate, merge, and write variables
		// each rank reads a subset of blocks
		for (size_t bi = rank; bi < nblocks; bi += nranks){
			const auto &block = blocks[bi];

			auto start = std::chrono::high_resolution_clock::now();
			if (verbose && rank == 0)
				std::cout << "step = " << step << ", block " << block.BlockID+1 << " out of " << nblocks << std::endl;


			//////////////////////////////////////////////////////////////////////////////////
			// read variables in my block

			// set block/subregion for the variables
			adios_connectivity.SetBlockSelection(block.BlockID);
			for (auto &coo : adios_coos) coo.SetBlockSelection(block.BlockID);
			for (auto &var : adios_vars) var.SetBlockSelection(block.BlockID);

			// allocate memory for reading variables
			std::vector<int64_t> ElementConnectivity;
			std::vector<std::vector<double>> coos(coo_names.size());
			std::vector<std::vector<double>> vars(var_names.size());

			// read variables
			bpReader.Get(adios_connectivity, ElementConnectivity, adios2::Mode::Sync);
			for (int i=0; i<coo_names.size(); i++) bpReader.Get(adios_coos[i], coos[i], adios2::Mode::Sync);
			for (int i=0; i<var_names.size(); i++) bpReader.Get(adios_vars[i], vars[i], adios2::Mode::Sync);
			bpReader.PerformGets();


			//////////////////////////////////////////////////////////////////////////////////
			// reindexing

			// Find index map: index of the node -> unique index of the node
			auto [ind_map, n_unique] = index_map_from_coos(coos);

			// reindex dofs to unique indices
			#ifdef _OPENMP
			#pragma omp parallel for schedule(static)
			#endif
			for(size_t i = 0; i < ElementConnectivity.size(); i++)
				ElementConnectivity[i] = ind_map[ElementConnectivity[i]];


			if (debug && rank == 0){
				// print Element Connectivity
				std::cout << "Element Connectivity:" << std::endl;
				int nnodes = n_nodes_in_element;
				for (int i=0; i<10*nnodes; i+=nnodes){
					for (int j=0; j<nnodes; j++) std::cout << " " << ElementConnectivity[i+j];
					std::cout << std::endl;
				}
				return 0;

				// // print variable at duplicate nodes
				// for (int vari=0; vari<vars.size(); vari++){
				// 	for (int i=0; i<inv_map[0].size(); i++){
				// 		std::cout << vars[vari][inv_map[0][i]] << " ";
				// 	}
				// 	std::cout << std::endl;
				// }
				// return 0;
			}

			//////////////////////////////////////////////////////////////////////////////////
			// merge variables, DG to CG

			std::vector<std::vector<double>> deduped_coos(coo_names.size());
			std::vector<std::vector<double>> deduped_vars(var_names.size());
			
			#ifdef _OPENMP
			#pragma omp parallel for schedule(static)
			#endif
			for (size_t i=0; i<coo_names.size(); i++) deduped_coos[i] = merge_values(ind_map, coos[i], n_unique);

			#ifdef _OPENMP
			#pragma omp parallel for schedule(static)
			#endif
			for (size_t i=0; i<var_names.size(); i++) deduped_vars[i] = merge_values(ind_map, vars[i], n_unique);


			//////////////////////////////////////////////////////////////////////////////////
			// save merged variables

			orig_adios_connectivity.SetSelection(adios2::Box<adios2::Dims>({}, {ElementConnectivity.size()}));
			bpWriter.Put<int64_t>(orig_adios_connectivity, ElementConnectivity.data(), adios2::Mode::Sync);
			
			for (int i=0; i<coo_names.size(); i++){
				out_adios_coos[i].SetSelection(adios2::Box<adios2::Dims>({}, {deduped_coos[i].size()}));
				bpWriter.Put<double>(out_adios_coos[i], deduped_coos[i].data(), adios2::Mode::Sync);
			}

			for (int i=0; i<var_names.size(); i++){
				out_adios_vars[i].SetSelection(adios2::Box<adios2::Dims>({}, {deduped_vars[i].size()}));
				bpWriter.Put<double>(out_adios_vars[i], deduped_vars[i].data(), adios2::Mode::Sync);
			}


			//////////////////////////////////////////////////////////////////////////////////
			// convert to CGNS standard: add element type and number of nodes in element

			if (convert_to_cgns_standard){	
				size_t num_elements = ElementConnectivity.size() / nodes_per_el;

				std::vector<int64_t> out_ElementConnectivity(ElementConnectivity.size());
				std::vector<int32_t> out_ElementNumNodes(num_elements);
				std::vector<uint8_t> out_ElementType(num_elements);

				size_t out_ElementConnectivity_size{0};
				for (size_t elem_id = 0, k = 0; k < ElementConnectivity.size(); k += nodes_per_el, ++elem_id) {
					auto* in_begin  = ElementConnectivity.data() + k;
					auto* in_end    = in_begin + nodes_per_el;
					auto* out_begin = out_ElementConnectivity.data() + out_ElementConnectivity_size;

					auto* out_end = std::unique_copy(in_begin, in_end, out_begin);
					int uniq_nodes_per_el = static_cast<int>(out_end - out_begin);

					out_ElementNumNodes[elem_id] = static_cast<int32_t>(uniq_nodes_per_el);
					out_ElementType[elem_id] = elem_type_from_num_nodes(uniq_nodes_per_el);
					out_ElementConnectivity_size += static_cast<size_t>(uniq_nodes_per_el);
				}

				out_adios_connectivity.SetSelection(adios2::Box<adios2::Dims>({}, {out_ElementConnectivity_size}));
				bpWriter.Put<int64_t>(out_adios_connectivity, out_ElementConnectivity.data(), adios2::Mode::Sync);

				out_adios_elemnumnodes.SetSelection(adios2::Box<adios2::Dims>({}, {out_ElementNumNodes.size()}));
				bpWriter.Put<int32_t>(out_adios_elemnumnodes, out_ElementNumNodes.data(), adios2::Mode::Sync);

				out_adios_elemtype.SetSelection(adios2::Box<adios2::Dims>({}, {out_ElementType.size()}));
				bpWriter.Put<uint8_t>(out_adios_elemtype, out_ElementType.data(), adios2::Mode::Sync);
			}

			bpWriter.PerformPuts();


			///////////////////////////////////////////////////////////////

			if (verbose && rank == 0){
				auto end = std::chrono::high_resolution_clock::now();
				auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

				total_nodes += coos[0].size();
				total_unique_nodes += n_unique;
				total_time += (double)duration.count() / 1e6;
				
				std::cout << "\tNum. orig./uniq. nodes = " << coos[0].size() << "/" << n_unique
						<< std::setprecision(3)
						<< ", redundancy = " <<  (double)coos[0].size()/(double)n_unique
						<< ", time = " << (double)duration.count() / 1e6 << " sec" << std::endl;
			}


			///////////////////////////////////////////////////////////////

			if (max_block_id>=0 and block.BlockID>=max_block_id) break;
		} // end for myBlockIdx

		bpReader.EndStep();	// end logical step
		bpWriter.EndStep();	// end logical step

	} // end for time step

	if (verbose && rank == 0){
		std::cout << "Rank " << rank << std::endl;
		std::cout << "Total num. orig. nodes = " << total_nodes << std::endl;
		std::cout << "Total num. uniq. nodes = " << total_unique_nodes << std::endl;
		std::cout << "Total redundancy = " << (double)total_nodes/(double)total_unique_nodes << std::endl;
		std::cout << "Total time = " << std::setprecision(3) << total_time << std::endl;
	}

	bpReader.Close();	// close engine
	bpWriter.Close();	// close engine


#if ADIOS2_USE_MPI
	MPI_Finalize();
#endif
	return 0;
}
