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


/**
 * Calculate global node connectivity
 */
template <typename Container>
auto ComputeGlobalConnectivity( const Container ElementConnectivity, size_t num_nodes, int nodes_per_el )
-> std::vector<std::set<size_t>>
{
	std::vector<std::set<size_t>> GlobalConnectivity(num_nodes);

	for(size_t k = 0; k < ElementConnectivity.size(); k+=nodes_per_el){
		for(int i = 0; i < nodes_per_el; i++){
			size_t n1 = ElementConnectivity[k+i];
			for(size_t j = 0; j < nodes_per_el; j++){
				size_t n2 = ElementConnectivity[k+j];
				if (n1!=n2){
					GlobalConnectivity[n1].insert(n2);
					GlobalConnectivity[n2].insert(n1);
				}
			}
		}
	}
	return GlobalConnectivity;
}


template <typename T>
auto reorder(const std::vector<T> &vec, const std::vector<int64_t> &order)
-> std::vector<T>
{
	std::vector<T> new_vec(vec.size());

	for (size_t i=0; i<vec.size(); i++)
		new_vec[i] = vec[order[i]];

	return new_vec;
}


/**
 * Calculate inverse index map: unique index of the node -> all indices of this node
 *
 * Inputs:
 * ------
 *  GlobalConnectivity - node connectivity: node id -> connected ids
 * 	coos - vector of x-y-z coordinate vectors
 */
template <typename Container>
auto find_node_order(const std::vector<std::set<size_t>> &GlobalConnectivity, const Container &coos)
-> std::vector<int64_t>
{
	size_t num_nodes = GlobalConnectivity.size();

	std::vector<bool>   node_notvisited(num_nodes, true);
	std::vector<int64_t> node_order(num_nodes, 0);

	// start with first node
	node_notvisited[0] = false;
	size_t curr_node = 0;
	size_t node_num_visited = 1;
	size_t notvisited_offset = 1;

	// walk through the graph using GlobalConnectivity
	// auto start = std::chrono::high_resolution_clock::now();
	while (node_num_visited<num_nodes){
		// find connected nodes that have not been visited yet
		std::vector<size_t> nns;
		for (const auto &n : GlobalConnectivity[curr_node])
			if (node_notvisited[n])
				nns.push_back(n);

		// if((node_num_visited-1)%1000==0){
		// 	std::cout << curr_node << ": ";
		// 	for (const auto & n : nns) std::cout << n << " ";
		// 	std::cout << std::endl;
		// }

		size_t next_node;
		if (nns.size()>0){
			double min_dst = 1.e15;
			for (const auto &candidate : nns){
				double dst = 0;
				for (auto &coo : coos)
					dst += std::pow(coo[curr_node]-coo[candidate],2);
				if (dst<min_dst){
					min_dst   = dst;
					next_node = candidate;
				}
			}
		}else{
			// find first not visited node
			next_node = std::find(node_notvisited.begin()+notvisited_offset, node_notvisited.end(), true) - node_notvisited.begin();
			notvisited_offset = next_node;
		}
		node_notvisited[next_node]   = false;
		node_order[node_num_visited] = next_node;

		curr_node = next_node;
		node_num_visited++;

		// auto end = std::chrono::high_resolution_clock::now();
		// auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
		// if(node_num_visited%1000==0){
		// 	start = std::chrono::high_resolution_clock::now();

		// 	std::cout << node_num_visited << " " << num_nodes << " " << curr_node << " " << (double)duration.count() / 1e6 << " " << (float)(node_num_visited)/num_nodes << std::endl;
		// }
	}

	return node_order;
}




///////////////////////////////////////////////////////////////////////////////


#endif