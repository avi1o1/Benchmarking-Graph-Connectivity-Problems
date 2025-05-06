#include <iostream>
#include <vector>
#include <fstream>
#include <chrono>
#include <string>
#include <utility>
#include <algorithm>
#include <unordered_map>

class Graph {
public:
    int num_nodes;
    std::vector<std::vector<int>> adj_list;
    
    // Constructor
    Graph(int n) : num_nodes(n), adj_list(n) {}
    
    // Add an edge to the graph
    void add_edge(int u, int v) {
        if (u >= 0 && u < num_nodes && v >= 0 && v < num_nodes && u != v) {
            // Check if edge already exists
            if (std::find(adj_list[u].begin(), adj_list[u].end(), v) == adj_list[u].end()) {
                adj_list[u].push_back(v);
                adj_list[v].push_back(u);
            }
        }
    }
    
    // Get node degree
    int degree(int node) const {
        return adj_list[node].size();
    }
};

// Union-Find (Disjoint Set) data structure for finding connected components
class UnionFind {
private:
    std::vector<int> parent;
    std::vector<int> rank;
    
public:
    // Constructor
    UnionFind(int n) : parent(n), rank(n, 0) {
        // Initialize each element as a separate set
        for (int i = 0; i < n; i++) {
            parent[i] = i;
        }
    }
    
    // Find with path compression
    int find(int x) {
        if (parent[x] != x) {
            parent[x] = find(parent[x]); // Path compression
        }
        return parent[x];
    }
    
    // Union by rank
    void unite(int x, int y) {
        int root_x = find(x);
        int root_y = find(y);
        
        if (root_x == root_y) return;
        
        // Union by rank
        if (rank[root_x] < rank[root_y]) {
            parent[root_x] = root_y;
        } else {
            parent[root_y] = root_x;
            if (rank[root_x] == rank[root_y]) {
                rank[root_x]++;
            }
        }
    }
    
    // Get all the components
    std::vector<std::vector<int>> get_components() {
        std::unordered_map<int, std::vector<int>> component_map;
        
        // Group elements by their root
        for (int i = 0; i < parent.size(); i++) {
            int root = find(i);
            component_map[root].push_back(i);
        }
        
        // Convert the map to a vector of components
        std::vector<std::vector<int>> components;
        for (const auto& pair : component_map) {
            components.push_back(pair.second);
        }
        
        return components;
    }
};

// Bridge finder - identifies bridges in a graph using DFS
class BridgeFinder {
private:
    Graph& graph;
    int time;
    std::vector<bool> visited;
    std::vector<int> discovery;
    std::vector<int> low;
    std::vector<std::pair<int, int>> bridges;
    
    // DFS to find bridges
    void dfs(int u, int parent) {
        visited[u] = true;
        discovery[u] = low[u] = ++time;
        
        for (int v : graph.adj_list[u]) {
            // If v is not visited
            if (!visited[v]) {
                dfs(v, u);
                
                // Check if the edge u-v is a bridge
                low[u] = std::min(low[u], low[v]);
                
                if (low[v] > discovery[u]) {
                    // This is a bridge
                    bridges.push_back(std::make_pair(u, v));
                }
            }
            // Update low value if v is already visited and not parent
            else if (v != parent) {
                low[u] = std::min(low[u], discovery[v]);
            }
        }
    }
    
public:
    // Constructor
    BridgeFinder(Graph& g) : graph(g), time(0), 
                            visited(g.num_nodes, false),
                            discovery(g.num_nodes, -1),
                            low(g.num_nodes, -1) {}
    
    // Find all bridges
    std::vector<std::pair<int, int>> find_bridges() {
        bridges.clear();
        
        // Call DFS for all unvisited vertices
        for (int i = 0; i < graph.num_nodes; i++) {
            if (!visited[i]) {
                dfs(i, -1);
            }
        }
        
        return bridges;
    }
};

// Load a graph from a file
Graph load_graph_from_file(const std::string& filename) {
    std::ifstream file(filename);
    if (!file) {
        std::cerr << "Error opening file " << filename << std::endl;
        exit(EXIT_FAILURE);
    }
    
    int num_nodes;
    file >> num_nodes;
    
    std::cout << "Initializing graph with " << num_nodes << " nodes..." << std::endl;
    Graph graph(num_nodes);
    
    int u, v;
    int edge_count = 0;
    std::cout << "Loading edges..." << std::endl;
    
    while (file >> u >> v) {
        // Convert 1-based input to 0-based indexing if needed
        graph.add_edge(u-1, v-1);  // Subtract 1 to convert to 0-based
        edge_count++;
        
        // Show progress for large files
        if (edge_count % 1000000 == 0) {
            std::cout << "Loaded " << (edge_count / 1000000) << " million edges..." << std::endl;
        }
    }
    
    std::cout << "Graph loaded: " << num_nodes << " nodes, " << edge_count << " edges" << std::endl;
    file.close();
    return graph;
}

// Print the biconnected components
void print_biconnected_components(const std::vector<std::vector<int>>& components, 
                                 const std::vector<std::pair<int, int>>& bridges) {
    std::cout << "Bridges (critical edges):" << std::endl;
    for (const auto& bridge : bridges) {
        std::cout << bridge.first << " -- " << bridge.second << std::endl;
    }
    std::cout << "Total bridges: " << bridges.size() << std::endl;
    
    std::cout << "\nBiconnected Components: " << components.size() << std::endl;
    
    // Sort components by size (largest first)
    std::vector<std::vector<int>> sorted_components = components;
    std::sort(sorted_components.begin(), sorted_components.end(), 
              [](const std::vector<int>& a, const std::vector<int>& b) {
                  return a.size() > b.size();
              });
    
    int max_to_print = std::min(20, static_cast<int>(sorted_components.size()));
    
    for (int i = 0; i < max_to_print; i++) {
        std::cout << "Component " << (i+1) << " (size " << sorted_components[i].size() << "): ";
        int max_nodes_to_print = std::min(10, static_cast<int>(sorted_components[i].size()));
        
        for (int j = 0; j < max_nodes_to_print; j++) {
            std::cout << sorted_components[i][j] << " ";
        }
        
        if (sorted_components[i].size() > max_nodes_to_print) {
            std::cout << "... (and " << (sorted_components[i].size() - max_nodes_to_print) << " more)";
        }
        std::cout << std::endl;
    }
    
    if (components.size() > max_to_print) {
        std::cout << "... (and " << (components.size() - max_to_print) << " more components)" << std::endl;
    }
}

int main(int argc, char *argv[]) {
    std::string filename;
    
    // Use command line argument for filename if provided
    if (argc > 1) {
        filename = argv[1];
    } else {
        filename = "datasets/custom.txt";
    }
    
    auto start_time = std::chrono::high_resolution_clock::now();
    
    std::cout << "Loading graph from " << filename << "..." << std::endl;
    Graph graph = load_graph_from_file(filename);
    
    std::cout << "Graph loaded from " << filename << " with " << graph.num_nodes << " nodes" << std::endl;
    std::cout << "Running sequential computation..." << std::endl;
    
    auto compute_start = std::chrono::high_resolution_clock::now();
    
    // Step 1: Find bridges
    BridgeFinder bridge_finder(graph);
    auto bridges = bridge_finder.find_bridges();
    
    // Step 2: Use Union-Find for biconnected components
    UnionFind uf(graph.num_nodes);
    
    // Create a map of bridges for quick lookup
    std::unordered_map<int, std::vector<int>> bridge_map;
    for (const auto& bridge : bridges) {
        bridge_map[bridge.first].push_back(bridge.second);
        bridge_map[bridge.second].push_back(bridge.first);
    }
    
    // Unite vertices connected by non-bridge edges
    for (int u = 0; u < graph.num_nodes; u++) {
        for (int v : graph.adj_list[u]) {
            if (u < v) { // Process each edge only once
                // Check if this edge is a bridge
                bool is_bridge = false;
                if (bridge_map.find(u) != bridge_map.end()) {
                    if (std::find(bridge_map[u].begin(), bridge_map[u].end(), v) != bridge_map[u].end()) {
                        is_bridge = true;
                    }
                }
                
                // If not a bridge, unite the vertices
                if (!is_bridge) {
                    uf.unite(u, v);
                }
            }
        }
    }
    
    // Get the biconnected components
    auto components = uf.get_components();
    
    auto compute_end = std::chrono::high_resolution_clock::now();
    auto compute_duration = std::chrono::duration_cast<std::chrono::microseconds>(compute_end - compute_start);
    
    std::cout << "Computation time: " << (compute_duration.count() / 1000000.0) << " seconds" << std::endl;
    print_biconnected_components(components, bridges);
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto total_duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
    
    std::cout << "Total execution time: " << (total_duration.count() / 1000000.0) << " seconds" << std::endl;
    
    return 0;
}

/**
 * Compilation and Execution Instructions:
 * 
 * To compile this program:
 * g++ -std=c++11 -o slota_madduri slota_madduri_seq.cpp
 * 
 * To run the program with a custom input file:
 * ./slota_madduri /path/to/your/graph/file.txt
 * 
 * To run with the default file path:
 * ./slota_madduri
 * 
 * Output:
 * The program outputs to the console (standard output):
 * 1. Information about loading the graph
 * 2. List of bridges (critical edges) found in the graph
 * 3. List of biconnected components found in the graph
 * 4. Timing information about the computation
 *
 * Note: Make sure the input file exists and has the correct format:
 * - First line: number of nodes
 * - Subsequent lines: edges as "u v" pairs
 */
