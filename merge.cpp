#include <adios2.h>
#include <vector>
#include <algorithm>    // std::find

#include <iostream>
#include <fstream>

#include <nlohmann/json.hpp>

using json = nlohmann::json;




/**
 * Calculate inverse index map: unique index of the node -> all indices of this node
 */
template <typename Container>
auto inverse_index_map( const Container &coos )
-> std::vector<std::vector<size_t>>
{
	int ndim = coos.size();

	// typedef std::vector<typename Container::value_type> node;
	typedef std::vector<double> node;
	typedef std::vector<size_t> indices;

	// collect all indices of each node
	size_t unique_ind = 0;
	node curr(ndim);
	std::map<node,indices> nodemap;
	for (size_t i = 0; i < coos[0].size(); i++){
		// ndim-dimensional node
		for (int d=0; d<ndim; d++) curr[d] = coos[d][i];

		// assign unique index
		if(nodemap[curr].empty()) nodemap[curr].push_back(unique_ind++);

		// original indices of the current unique index
		nodemap[curr].push_back(i);
	}

	// map unique indices of the nodes to all indices of this node
	std::vector<indices> indmap(nodemap.size());
	for (const auto &nm : nodemap){
		unique_ind = nm.second[0];
		indmap[unique_ind] = indices{nm.second.begin()+1,nm.second.end()};
	}

	return indmap;
}


/**
 * Convert inverse index map to the forward map: index of the node -> unique index of the node
 */
auto index_map_from_inverse_map(std::vector<std::vector<size_t>> &inv_map, size_t n)
-> std::vector<size_t>
{
	std::vector<size_t> forward_map(n);

	for (size_t i=0; i<inv_map.size(); i++)
		for (const auto &j : inv_map[i])
			forward_map[j] = i;

	return forward_map;
}


/**
 * Convert inverse map to serial vector by using the format:
 * num_undices, unique_index, indices,
 * num_undices, unique_index, indices,
 * ...
 */
auto serialize_inv_map(std::map<size_t, std::vector<size_t>> &inv_map, size_t n)
-> std::vector<size_t>
{
	std::vector<size_t> serial_inv_map(inv_map.size()+n);

	size_t i = 0;
	for (const auto &si : inv_map){
		serial_inv_map[i++] = si.second.size();
		for (const auto &v : si.second){
			serial_inv_map[i++] = v;
		}
	}

	return serial_inv_map;
}



template <typename T>
auto merge_values(std::vector<std::vector<size_t>> &inv_map, std::vector<T> &values)
-> std::vector<T>
{
	std::vector<T> merged_values(inv_map.size());

	for (size_t i=0; i<inv_map.size(); i++)
	{
		T merged_val = 0;
		for (const auto & n : inv_map[i]) merged_val += values[n];
		merged_values[i] = merged_val / inv_map[i].size();
	}

	return merged_values;
}



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


	// // ADIOS compression
	// adios2::Operator op = adios.DefineOperator("mgardplus", "mgardplus");


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
	// adios2::Variable<size_t> inv_map_out = writer_io.DefineVariable<size_t>("inv_map", {}, {}, {adios2::UnknownDim});

	adios2::Variable<int64_t> out_adios_connectivity = writer_io.DefineVariable<int64_t>(connectivity_name, {}, {}, {adios2::UnknownDim});

	std::vector<adios2::Variable<double>> out_adios_coos(coo_names.size());
	for (int i=0; i<coo_names.size(); i++)
		out_adios_coos[i] = writer_io.DefineVariable<double>(coo_names[i], {}, {}, {adios2::UnknownDim});

	std::vector<adios2::Variable<double>> out_adios_vars(var_names.size());
	for (int i=0; i<var_names.size(); i++)
		out_adios_vars[i] = writer_io.DefineVariable<double>(var_names[i], {}, {}, {adios2::UnknownDim});


	///////////////////////////////////////////////////////////////////////////


	// number of blocks/subdomains
	auto blocks = bpReader.BlocksInfo(adios_connectivity, 0);
	size_t nblocks = blocks.size();

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
		// reindexing

		// Calculate inverse index map: unique linear index of the node -> all indices of this (duplicate) node
		auto inv_map = inverse_index_map(coos);

		// Convert inverse index map to the forward map: index of the node -> unique index of the node
		auto ind_map = index_map_from_inverse_map(inv_map, coos[0].size());

		// // Convert inverse map to serial vector for saving in .bp file
		// auto serial_inv_map = serialize_inv_map(inv_map, coos[0].size());

		// reindex dofs to unique indices
		for(size_t i = 0; i < ElementConnectivity.size(); i++)
			ElementConnectivity[i] = ind_map[ElementConnectivity[i]];
		// // print Element Connectivity
		// {
		// std::cout << "\tElement Connectivity:" << std::endl;
		// int nnodes = config["num_nodes_in_bp_element"];
		// for (int i=0; i<10*nnodes; i+=nnodes)
		// {
		// 	for (int j=0; j<nnodes; j++) std::cout << " " << ElementConnectivity[i+j];
		// 	std::cout << std::endl;
		// }
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

		std::cout << "\tNum. orig. points = " << coos[0].size() << std::endl;
		std::cout << "\tNum. uniq. points = " << inv_map.size() << std::endl;

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
	bpReader.Close();	// close engine
	bpWriter.Close();	// close engine
}
