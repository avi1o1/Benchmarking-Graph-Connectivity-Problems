/*
 * To compile this file with OpenMP support, use:
 * gcc -O3 -fopenmp tarjan_vishkin_par.c -o tarjan_vishkin_par && ./tarjan_vishkin_par
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>
#include <omp.h>

// Adjacency list node
typedef struct AdjListNode {
    int vertex;
    int edge_id;  // To track edge IDs for the algorithm
    struct AdjListNode* next;
} AdjListNode;

// Graph structure with adjacency list
typedef struct {
    int num_nodes;
    int num_edges;
    AdjListNode** adjacency_list;
} Graph;

// Structure to represent an edge
typedef struct {
    int src, dest;
} Edge;

// Structure to represent a biconnected component
typedef struct {
    int* edges;
    int size;
    int capacity;
} BiconnectedComponent;

// Structure to hold multiple biconnected components
typedef struct {
    BiconnectedComponent* components;
    int count;
    int capacity;
    Edge* edges; // Store edges for later printing
} ComponentList;

// Structure for UNION-FIND operations
typedef struct {
    int* parent;
    int* rank;
    int size;
} DisjointSet;

// Initialize disjoint set
void init_disjoint_set(DisjointSet* ds, int size) {
    ds->size = size;
    ds->parent = (int*)malloc(size * sizeof(int));
    ds->rank = (int*)calloc(size, sizeof(int));
    
    if (!ds->parent || !ds->rank) {
        fprintf(stderr, "Memory allocation failed for disjoint set\n");
        exit(EXIT_FAILURE);
    }
    
    #pragma omp parallel for
    for (int i = 0; i < size; i++) {
        ds->parent[i] = i;  // Each element is its own parent initially
    }
}

// Find root of the set with path compression
int find_parallel(DisjointSet* ds, int x) {
    if (ds->parent[x] != x) {
        int root = find_parallel(ds, ds->parent[x]);
        ds->parent[x] = root;
        return root;
    }
    return ds->parent[x];
}

// Find root of the set for parallel processing (no recursive calls)
int find(DisjointSet* ds, int x) {
    int root = x;
    
    // Find the root
    while (ds->parent[root] != root) {
        root = ds->parent[root];
    }
    
    // Compress path
    while (ds->parent[x] != root) {
        int next = ds->parent[x];
        ds->parent[x] = root;
        x = next;
    }
    
    return root;
}

// Union two sets with atomic operations for thread safety
void union_sets_atomic(DisjointSet* ds, int x, int y) {
    int root_x, root_y;
    
    do {
        root_x = find(ds, x);
        root_y = find(ds, y);
        
        if (root_x == root_y) return;
        
        // Union by rank with atomic compare and swap
        if (ds->rank[root_x] < ds->rank[root_y]) {
            // Atomically update parent
            #pragma omp atomic write
            ds->parent[root_x] = root_y;
        } else if (ds->rank[root_x] > ds->rank[root_y]) {
            // Atomically update parent
            #pragma omp atomic write
            ds->parent[root_y] = root_x;
        } else {
            // Equal ranks, update parent and increment rank
            #pragma omp atomic write
            ds->parent[root_y] = root_x;
            
            #pragma omp atomic update
            ds->rank[root_x]++;
        }
    } while (ds->parent[root_x] != root_y && ds->parent[root_y] != root_x);
}

// Thread-local union operation
void union_sets(DisjointSet* ds, int x, int y) {
    int root_x = find(ds, x);
    int root_y = find(ds, y);
    
    if (root_x == root_y) return;
    
    // Union by rank
    if (ds->rank[root_x] < ds->rank[root_y]) {
        ds->parent[root_x] = root_y;
    } else if (ds->rank[root_x] > ds->rank[root_y]) {
        ds->parent[root_y] = root_x;
    } else {
        ds->parent[root_y] = root_x;
        ds->rank[root_x]++;
    }
}

// Free disjoint set
void free_disjoint_set(DisjointSet* ds) {
    free(ds->parent);
    free(ds->rank);
    ds->size = 0;
}

// Initialize an empty graph with n nodes
void init_graph(Graph* graph, int n) {
    graph->num_nodes = n;
    graph->num_edges = 0;
    graph->adjacency_list = (AdjListNode**)malloc(n * sizeof(AdjListNode*));
    if (!graph->adjacency_list) {
        fprintf(stderr, "Failed to allocate adjacency list\n");
        exit(EXIT_FAILURE);
    }
    
    #pragma omp parallel for
    for (int i = 0; i < n; i++) {
        graph->adjacency_list[i] = NULL;
    }
}

// Create a new adjacency list node
AdjListNode* new_adj_list_node(int dest, int edge_id) {
    AdjListNode* new_node = (AdjListNode*)malloc(sizeof(AdjListNode));
    if (!new_node) {
        fprintf(stderr, "Failed to allocate adjacency list node\n");
        exit(EXIT_FAILURE);
    }
    new_node->vertex = dest;
    new_node->edge_id = edge_id;
    new_node->next = NULL;
    return new_node;
}

// Add an edge to the graph (thread-safe)
void add_edge(Graph *graph, int u, int v) {
    if (u >= 0 && u < graph->num_nodes && v >= 0 && v < graph->num_nodes && u != v) {
        // Using a critical section to ensure thread safety
        #pragma omp critical
        {
            // Check if edge already exists
            AdjListNode* temp = graph->adjacency_list[u];
            bool exists = false;
            
            while (temp) {
                if (temp->vertex == v) {
                    exists = true;
                    break;
                }
                temp = temp->next;
            }
            
            if (!exists) {
                int edge_id = graph->num_edges;
                
                // Add edge u -> v
                AdjListNode* new_node = new_adj_list_node(v, edge_id);
                new_node->next = graph->adjacency_list[u];
                graph->adjacency_list[u] = new_node;
                
                // Add edge v -> u
                new_node = new_adj_list_node(u, edge_id);
                new_node->next = graph->adjacency_list[v];
                graph->adjacency_list[v] = new_node;
                
                graph->num_edges++;
            }
        }
    }
}

// Initialize a biconnected component
void init_component(BiconnectedComponent* component, int initial_capacity) {
    component->size = 0;
    component->capacity = initial_capacity;
    component->edges = (int*)malloc(initial_capacity * sizeof(int));
    if (!component->edges) {
        fprintf(stderr, "Failed to allocate component edges\n");
        exit(EXIT_FAILURE);
    }
}

// Add an edge to a biconnected component
void add_edge_to_component(BiconnectedComponent* component, int edge_id) {
    if (component->size == component->capacity) {
        component->capacity *= 2;
        component->edges = (int*)realloc(component->edges, component->capacity * sizeof(int));
        if (!component->edges) {
            fprintf(stderr, "Failed to reallocate component edges\n");
            exit(EXIT_FAILURE);
        }
    }
    component->edges[component->size++] = edge_id;
}

// Initialize component list
void init_component_list(ComponentList* list, int initial_capacity) {
    list->count = 0;
    list->capacity = initial_capacity;
    list->components = (BiconnectedComponent*)malloc(initial_capacity * sizeof(BiconnectedComponent));
    list->edges = NULL; // Initialize the edges pointer
    if (!list->components) {
        fprintf(stderr, "Failed to allocate component list\n");
        exit(EXIT_FAILURE);
    }
}

// Add a component to the component list
void add_component_to_list(ComponentList* list, int* edges, int size) {
    // Critical section to ensure thread safety when adding components
    #pragma omp critical
    {
        if (list->count == list->capacity) {
            list->capacity *= 2;
            list->components = (BiconnectedComponent*)realloc(list->components, 
                                                             list->capacity * sizeof(BiconnectedComponent));
            if (!list->components) {
                fprintf(stderr, "Failed to reallocate component list\n");
                exit(EXIT_FAILURE);
            }
        }
        
        init_component(&list->components[list->count], size);
        for (int i = 0; i < size; i++) {
            add_edge_to_component(&list->components[list->count], edges[i]);
        }
        list->count++;
    }
}

// Free component list memory
void free_component_list(ComponentList* list) {
    for (int i = 0; i < list->count; i++) {
        free(list->components[i].edges);
    }
    free(list->components);
    if (list->edges) {
        free(list->edges);
    }
    list->edges = NULL;
    list->components = NULL;
    list->count = 0;
    list->capacity = 0;
}

// DFS for computing preorder and low values
void dfs(Graph* graph, int u, int parent, int* preorder, int* low, Edge* edges, int* visited, int* time) {
    visited[u] = 1;
    
    // Atomic update for the time variable
    int current_time;
    #pragma omp atomic capture
    {
        current_time = *time;
        (*time)++;
    }
    
    preorder[u] = low[u] = current_time;
    
    AdjListNode* temp = graph->adjacency_list[u];
    
    while (temp) {
        int v = temp->vertex;
        
        if (v != parent) {  // Avoid parent edge
            if (!visited[v]) {
                // Store the edge information
                edges[temp->edge_id].src = u;
                edges[temp->edge_id].dest = v;
                
                // Recursive DFS
                dfs(graph, v, u, preorder, low, edges, visited, time);
                
                // Update low value
                if (low[v] < low[u]) {
                    low[u] = low[v];
                }
                
                // Check if this is a back edge
                if (low[v] >= preorder[u]) {
                    // u is an articulation point, except for the root
                    // We don't need to track articulation points here
                }
            } else if (preorder[v] < preorder[u]) {
                // This is a back edge
                edges[temp->edge_id].src = u;
                edges[temp->edge_id].dest = v;
                
                // Update low value
                if (preorder[v] < low[u]) {
                    low[u] = preorder[v];
                }
            }
        }
        
        temp = temp->next;
    }
}

// Parallel DFS starter for multiple connected components
void parallel_dfs_starter(Graph* graph, int* preorder, int* low, Edge* edges, int* visited, int* time) {
    int n = graph->num_nodes;
    bool* processed = (bool*)calloc(n, sizeof(bool));
    
    if (!processed) {
        fprintf(stderr, "Memory allocation failed for processed array\n");
        exit(EXIT_FAILURE);
    }
    
    // Process each node in parallel
    #pragma omp parallel
    {
        // Each thread tries to start DFS from unvisited nodes
        while (true) {
            int start_node = -1;
            
            // Find an unprocessed node
            #pragma omp critical
            {
                for (int i = 0; i < n; i++) {
                    if (!processed[i] && !visited[i]) {
                        start_node = i;
                        processed[i] = true;
                        break;
                    }
                }
            }
            
            // If no more nodes to process, exit the loop
            if (start_node == -1) {
                break;
            }
            
            // Process this connected component using DFS
            dfs(graph, start_node, -1, preorder, low, edges, visited, time);
        }
    }
    
    free(processed);
}

// Process edges for biconnected components in parallel
void process_edges_parallel(Graph* graph, int* preorder, int* low, Edge* edges, DisjointSet* ds) {
    int n = graph->num_nodes;
    
    // Process all edges in parallel
    #pragma omp parallel for schedule(dynamic, 128)
    for (int u = 0; u < n; u++) {
        AdjListNode* temp = graph->adjacency_list[u];
        
        while (temp) {
            int v = temp->vertex;
            int edge_id = temp->edge_id;
            
            if (edges[edge_id].src == u) { // Process each edge only once
                if (low[v] >= preorder[u] && preorder[u] < preorder[v]) {
                    // This is a bridge or part of a cut vertex, don't merge
                } else {
                    // Find other edges that should be in the same biconnected component
                    AdjListNode* other = graph->adjacency_list[u];
                    while (other) {
                        if (other->edge_id != edge_id) {
                            int w = other->vertex;
                            if (edges[other->edge_id].src == u) {
                                if ((preorder[w] < preorder[u] && preorder[w] < preorder[v] && low[v] <= preorder[w]) ||
                                    (preorder[v] < preorder[u] && preorder[v] < preorder[w] && low[w] <= preorder[v])) {
                                    // These edges are in the same biconnected component
                                    union_sets_atomic(ds, edge_id, other->edge_id);
                                }
                            }
                        }
                        other = other->next;
                    }
                }
            }
            
            temp = temp->next;
        }
    }
}

// Tarjan-Vishkin algorithm for finding biconnected components - parallel version
void tarjan_vishkin_biconnected_components_parallel(Graph* graph, ComponentList* result) {
    int n = graph->num_nodes;
    int m = graph->num_edges;
    
    // Arrays for DFS traversal
    int* preorder = (int*)malloc(n * sizeof(int));
    int* low = (int*)malloc(n * sizeof(int));
    int* visited = (int*)calloc(n, sizeof(int));
    Edge* edges = (Edge*)calloc(m, sizeof(Edge)); // Changed to calloc to ensure initialization
    
    if (!preorder || !low || !visited || !edges) {
        fprintf(stderr, "Memory allocation failed for DFS arrays\n");
        exit(EXIT_FAILURE);
    }
    
    // Initialize arrays in parallel
    #pragma omp parallel for
    for (int i = 0; i < n; i++) {
        preorder[i] = 0;
        low[i] = 0;
        visited[i] = 0;
    }
    
    // Initialize DFS time
    int time = 0;
    
    // Perform DFS for each unvisited node with parallel starters
    parallel_dfs_starter(graph, preorder, low, edges, visited, &time);
    
    // Use disjoint-set to find biconnected components
    DisjointSet ds;
    init_disjoint_set(&ds, m);
    
    // Process edges based on low values in parallel
    process_edges_parallel(graph, preorder, low, edges, &ds);
    
    // Count the number of biconnected components
    int* component_id = (int*)malloc(m * sizeof(int));
    int component_count = 0;
    
    #pragma omp parallel
    {
        // First, have all threads find the root for each edge in parallel
        #pragma omp for
        for (int i = 0; i < m; i++) {
            component_id[i] = find(&ds, i);
        }
        
        // One thread counts the components
        #pragma omp single
        {
            int* is_root = (int*)calloc(m, sizeof(int));
            for (int i = 0; i < m; i++) {
                if (component_id[i] == i) {
                    is_root[i] = 1;
                }
            }
            
            for (int i = 0; i < m; i++) {
                if (is_root[i]) {
                    component_count++;
                }
            }
            
            free(is_root);
        }
    }
    
    // Create arrays for each component
    int** component_edges = (int**)malloc(component_count * sizeof(int*));
    int* component_sizes = (int*)calloc(component_count, sizeof(int));
    int* component_capacities = (int*)malloc(component_count * sizeof(int));
    int* component_map = (int*)malloc(m * sizeof(int));
    
    if (!component_edges || !component_sizes || !component_capacities || !component_map) {
        fprintf(stderr, "Memory allocation failed for component arrays\n");
        exit(EXIT_FAILURE);
    }
    
    // Initialize component maps
    int idx = 0;
    for (int i = 0; i < m; i++) {
        if (component_id[i] == i) {
            component_map[i] = idx++;
            component_capacities[component_map[i]] = 10;  // Initial capacity
            component_edges[component_map[i]] = (int*)malloc(10 * sizeof(int));
            if (!component_edges[component_map[i]]) {
                fprintf(stderr, "Memory allocation failed for component edges\n");
                exit(EXIT_FAILURE);
            }
        }
    }
    
    // Group edges by biconnected component
    #pragma omp parallel
    {
        // Thread-private counter arrays
        int* local_sizes = (int*)calloc(component_count, sizeof(int));
        
        // Each thread counts its contribution to each component
        #pragma omp for schedule(static)
        for (int i = 0; i < m; i++) {
            int root = component_id[i];
            int component = component_map[root];
            local_sizes[component]++;
        }
        
        // Allocate space in the global arrays
        #pragma omp critical
        {
            for (int c = 0; c < component_count; c++) {
                if (component_sizes[c] + local_sizes[c] > component_capacities[c]) {
                    while (component_sizes[c] + local_sizes[c] > component_capacities[c]) {
                        component_capacities[c] *= 2;
                    }
                    component_edges[c] = (int*)realloc(component_edges[c], 
                                                    component_capacities[c] * sizeof(int));
                    if (!component_edges[c]) {
                        fprintf(stderr, "Memory allocation failed for component edges\n");
                        exit(EXIT_FAILURE);
                    }
                }
            }
        }
        
        // Add edges to components using thread-local counters
        #pragma omp barrier
        
        int* local_offsets = (int*)malloc(component_count * sizeof(int));
        if (!local_offsets) {
            fprintf(stderr, "Memory allocation failed for local offsets\n");
            exit(EXIT_FAILURE);
        }
        
        #pragma omp critical
        {
            for (int c = 0; c < component_count; c++) {
                local_offsets[c] = component_sizes[c];
                component_sizes[c] += local_sizes[c];
            }
        }
        
        // Now each thread has its own portion of the arrays to fill
        #pragma omp for schedule(static)
        for (int i = 0; i < m; i++) {
            int root = component_id[i];
            int component = component_map[root];
            int offset = local_offsets[component]++;
            if (offset < component_capacities[component]) {
                component_edges[component][offset] = i;
            }
        }
        
        free(local_sizes);
        free(local_offsets);
    }
    
    // Add components to result
    init_component_list(result, component_count);
    
    for (int i = 0; i < component_count; i++) {
        add_component_to_list(result, component_edges[i], component_sizes[i]);
    }
    
    // Clean up
    for (int i = 0; i < component_count; i++) {
        free(component_edges[i]);
    }
    
    free(preorder);
    free(low);
    free(visited);
    // Store the edges for later printing
    result->edges = edges; // Pass ownership of edges to result
    
    free(component_id);
    free(component_edges);
    free(component_sizes);
    free(component_capacities);
    free(component_map);
    free_disjoint_set(&ds);
}

// Print biconnected components in a human-readable format
void print_biconnected_components(ComponentList* components, Edge* edges) {
    printf("Found %d biconnected components:\n", components->count);
    
    // Find max node id to allocate properly sized array
    int max_node_id = 0;
    for (int i = 0; i < components->count; i++) {
        for (int j = 0; j < components->components[i].size; j++) {
            int edge_id = components->components[i].edges[j];
            if (edge_id >= 0 && edges[edge_id].src > max_node_id) {
                max_node_id = edges[edge_id].src;
            }
            if (edge_id >= 0 && edges[edge_id].dest > max_node_id) {
                max_node_id = edges[edge_id].dest;
            }
        }
    }
    
    // Sort components by size (largest first)
    for (int i = 0; i < components->count - 1; i++) {
        for (int j = i + 1; j < components->count; j++) {
            if (components->components[i].size < components->components[j].size) {
                BiconnectedComponent temp = components->components[i];
                components->components[i] = components->components[j];
                components->components[j] = temp;
            }
        }
    }
    
    // Print each component in a standardized format matching slota_madduri
    for (int i = 0; i < components->count; i++) {
        printf("Component %d: Size: %d, Nodes: ", i + 1, components->components[i].size);
        
        // Create an array to track which nodes are in this component - dynamically sized
        bool* node_in_component = (bool*)calloc(max_node_id + 1, sizeof(bool));
        if (!node_in_component) {
            fprintf(stderr, "Memory allocation failed for node tracking\n");
            return;
        }
        
        // Mark all nodes that appear in this component
        for (int j = 0; j < components->components[i].size; j++) {
            int edge_id = components->components[i].edges[j];
            // Ensure edge_id is valid
            if (edge_id >= 0) {
                int src = edges[edge_id].src;
                int dest = edges[edge_id].dest;
                if (src >= 0 && src <= max_node_id) {
                    node_in_component[src] = true;
                }
                if (dest >= 0 && dest <= max_node_id) {
                    node_in_component[dest] = true;
                }
            }
        }
        
        // Print nodes in order
        for (int node = 0; node <= max_node_id; node++) {
            if (node_in_component[node]) {
                printf("%d ", node);
            }
        }
        printf("\n");
        
        free(node_in_component);
    }
}

// Free graph memory
void free_graph(Graph* graph) {
    for (int i = 0; i < graph->num_nodes; i++) {
        AdjListNode* current = graph->adjacency_list[i];
        while (current) {
            AdjListNode* temp = current;
            current = current->next;
            free(temp);
        }
    }
    
    free(graph->adjacency_list);
    free(graph);
}

// Load a graph from a file
Graph* load_graph_from_file(const char* filename) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        fprintf(stderr, "Error opening file %s\n", filename);
        return NULL;
    }
    
    Graph* graph = (Graph*)malloc(sizeof(Graph));
    if (!graph) {
        fprintf(stderr, "Memory allocation failed for graph\n");
        fclose(file);
        return NULL;
    }
    
    int num_nodes;
    if (fscanf(file, "%d", &num_nodes) != 1) {
        fprintf(stderr, "Error reading number of nodes\n");
        free(graph);
        fclose(file);
        return NULL;
    }
    
    printf("Initializing graph with %d nodes...\n", num_nodes);
    init_graph(graph, num_nodes);
    
    char line[256];
    // Check return value of fgets to avoid warning
    if (fgets(line, sizeof(line), file) == NULL) {
        fprintf(stderr, "Error reading file header\n");
        free_graph(graph);
        fclose(file);
        return NULL;
    }
    
    int edge_count = 0;
    printf("Loading edges...\n");
    
    while (fgets(line, sizeof(line), file)) {
        int u, v;
        if (sscanf(line, "%d %d", &u, &v) == 2) {
            // Convert 1-based input to 0-based indexing
            add_edge(graph, u-1, v-1);  // Subtract 1 to convert to 0-based
            edge_count++;
            
            // Show progress for large files
            if (edge_count % 1000000 == 0) {
                printf("Loaded %d million edges...\n", edge_count / 1000000);
            }
        }
    }
    
    printf("Graph loaded: %d nodes, %d edges\n", num_nodes, edge_count);
    fclose(file);
    return graph;
}

int main(int argc, char *argv[]) {
    const char* filename;
    
    // Use command line argument for filename if provided
    if (argc > 1) {
        filename = argv[1];
    } else {
        filename = "../datasets/small.txt";
    }
    
    clock_t start_time = clock();
    
    printf("Loading graph from %s...\n", filename);
    Graph* graph = load_graph_from_file(filename);
    
    if (!graph) {
        fprintf(stderr, "Failed to load graph from %s\n", filename);
        return 1;
    }
    
    printf("Graph loaded from %s with %d nodes\n", filename, graph->num_nodes);
    
    ComponentList components;
    printf("Running parallel Tarjan-Vishkin computation with %d threads...\n", omp_get_max_threads());
    
    clock_t compute_start = clock();
    tarjan_vishkin_biconnected_components_parallel(graph, &components);
    clock_t compute_end = clock();
    
    printf("Computation time: %.5f seconds\n", 
           (double)(compute_end - compute_start) / CLOCKS_PER_SEC);
    
    print_biconnected_components(&components, components.edges);
    
    // Clean up
    free_component_list(&components);
    free_graph(graph);
    
    clock_t end_time = clock();
    printf("Total execution time: %.5f seconds\n", 
           (double)(end_time - start_time) / CLOCKS_PER_SEC);
    
    return 0;
}