#ifndef NONUNIFORMMAP
#define NONUNIFORMMAP

#include <vector>
#include <map>
#include <algorithm>    // std::find




///////////////////////////////////////////////////////////////////////////////
// merging

/**
 * Calculate inverse index map: unique index of the node -> all indices of this node
 *
 * Inputs:
 * ------
 * 	coos - vector of x-y-z coordinate vectors
 */
template <typename Container>
std::vector<std::vector<size_t>> inverse_index_map( const Container &coos )
{
	int ndim = coos.size();

	typedef std::vector<typename Container::value_type::value_type> node;
	// typedef std::vector<double> node;
	typedef std::vector<size_t> indices;

	// collect all indices of each unique node
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
std::vector<size_t> index_map_from_inverse_map(std::vector<std::vector<size_t>> &inv_map, size_t n)
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
std::vector<size_t> serialize_inv_map(std::map<size_t, std::vector<size_t>> &inv_map, size_t n)
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
std::vector<T> merge_values(std::vector<std::vector<size_t>> &inv_map, std::vector<T> &values)
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


// auto adjacency_matrix(std::vector<int64_t> &connectivity, int nnodes)
// -> std::map<int64_t, std::set<size_t>>
// {
// 	std::map<int64_t, std::vector<int64_t>> adj_mat;

// 	// for (size_t el=0; el<connectivity.size()/nnodes; el++){
// 	for (size_t i=0; i<connectivity.size(); i+=nnodes){
// 		for (size_t j=i; j<nnodes; j++){
// 			adj_mat[i].insert(j);
// 			adj_mat[j].insert(i);
// 		}
// 	}

// 	return adj_mat;
// }




///////////////////////////////////////////////////////////////////////////////


#endif