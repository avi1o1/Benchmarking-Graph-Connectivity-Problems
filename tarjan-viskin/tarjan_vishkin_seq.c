/**
 * Implementation of the Tarjan-Vishkin algorithm for finding biconnected components in graphs.
 * This is a sequential C implementation.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>

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
    
    for (int i = 0; i < size; i++) {
        ds->parent[i] = i;  // Each element is its own parent initially
    }
}

// Find root of the set with path compression
int find(DisjointSet* ds, int x) {
    if (ds->parent[x] != x) {
        ds->parent[x] = find(ds, ds->parent[x]);  // Path compression
    }
    return ds->parent[x];
}

// Union two sets by rank
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

// Add an edge to the graph
void add_edge(Graph *graph, int u, int v) {
    if (u >= 0 && u < graph->num_nodes && v >= 0 && v < graph->num_nodes && u != v) {
        // Check if edge already exists
        AdjListNode* temp = graph->adjacency_list[u];
        while (temp) {
            if (temp->vertex == v) return; // Edge already exists
            temp = temp->next;
        }
        
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
    if (!list->components) {
        fprintf(stderr, "Failed to allocate component list\n");
        exit(EXIT_FAILURE);
    }
}

// Add a component to the component list
void add_component_to_list(ComponentList* list, int* edges, int size) {
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

// Free component list memory
void free_component_list(ComponentList* list) {
    for (int i = 0; i < list->count; i++) {
        free(list->components[i].edges);
    }
    free(list->components);
    list->components = NULL;
    list->count = 0;
    list->capacity = 0;
}

// DFS for computing preorder and low values
void dfs(Graph* graph, int u, int parent, int* preorder, int* low, Edge* edges, int* visited, int* time) {
    visited[u] = 1;
    preorder[u] = low[u] = (*time)++;
    
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

// Tarjan-Vishkin algorithm for finding biconnected components
void tarjan_vishkin_biconnected_components(Graph* graph, ComponentList* result) {
    int n = graph->num_nodes;
    int m = graph->num_edges;
    
    // Arrays for DFS traversal
    int* preorder = (int*)malloc(n * sizeof(int));
    int* low = (int*)malloc(n * sizeof(int));
    int* visited = (int*)calloc(n, sizeof(int));
    Edge* edges = (Edge*)malloc(m * sizeof(Edge));
    
    if (!preorder || !low || !visited || !edges) {
        fprintf(stderr, "Memory allocation failed for DFS arrays\n");
        exit(EXIT_FAILURE);
    }
    
    // Initialize DFS time
    int time = 0;
    
    // Perform DFS for each unvisited node
    for (int i = 0; i < n; i++) {
        if (!visited[i]) {
            dfs(graph, i, -1, preorder, low, edges, visited, &time);
        }
    }
    
    // Use disjoint-set to find biconnected components
    DisjointSet ds;
    init_disjoint_set(&ds, m);
    
    // Process edges based on low values
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
                                    union_sets(&ds, edge_id, other->edge_id);
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
    
    // Count the number of biconnected components
    int* component_id = (int*)malloc(m * sizeof(int));
    int component_count = 0;
    
    for (int i = 0; i < m; i++) {
        int root = find(&ds, i);
        component_id[i] = root;
        if (root == i) {
            component_count++;
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
        if (find(&ds, i) == i) {
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
    for (int i = 0; i < m; i++) {
        int root = find(&ds, i);
        int component = component_map[root];
        
        if (component_sizes[component] == component_capacities[component]) {
            component_capacities[component] *= 2;
            component_edges[component] = (int*)realloc(component_edges[component], 
                                                    component_capacities[component] * sizeof(int));
            if (!component_edges[component]) {
                fprintf(stderr, "Memory allocation failed for component edges\n");
                exit(EXIT_FAILURE);
            }
        }
        
        component_edges[component][component_sizes[component]++] = i;
    }
    
    // Add components to result
    init_component_list(result, component_count);
    for (int i = 0; i < component_count; i++) {
        add_component_to_list(result, component_edges[i], component_sizes[i]);
        free(component_edges[i]);
    }
    
    // Clean up
    free(preorder);
    free(low);
    free(visited);
    free(edges);
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
    
    // Print each component (limit output for very large results)
    int max_to_print = (components->count > 100) ? 100 : components->count;
    for (int i = 0; i < max_to_print; i++) {
        printf("Component %d: Size: %d, Edges: ", i + 1, components->components[i].size);
        // Print only a few edges for large components
        int edges_to_print = (components->components[i].size > 10) ? 10 : components->components[i].size;
        for (int j = 0; j < edges_to_print; j++) {
            int edge_id = components->components[i].edges[j];
            printf("(%d-%d) ", edges[edge_id].src, edges[edge_id].dest);
        }
        if (edges_to_print < components->components[i].size) {
            printf("... (and %d more)", components->components[i].size - edges_to_print);
        }
        printf("\n");
    }
    
    if (max_to_print < components->count) {
        printf("... (showing only %d out of %d components)\n", max_to_print, components->count);
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
    fgets(line, sizeof(line), file);  // Skip the first line (number of nodes)
    
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
    
    // Prepare edge array for the results
    Edge* edges = (Edge*)malloc(graph->num_edges * sizeof(Edge));
    if (!edges) {
        fprintf(stderr, "Memory allocation failed for edges array\n");
        free_graph(graph);
        return 1;
    }
    
    ComponentList components;
    printf("Running sequential Tarjan-Vishkin computation...\n");
    
    clock_t compute_start = clock();
    tarjan_vishkin_biconnected_components(graph, &components);
    clock_t compute_end = clock();
    
    printf("Computation time: %.5f seconds\n", 
           (double)(compute_end - compute_start) / CLOCKS_PER_SEC);
    
    print_biconnected_components(&components, edges);
    
    // Clean up
    free_component_list(&components);
    free(edges);
    free_graph(graph);
    
    clock_t end_time = clock();
    printf("Total execution time: %.5f seconds\n", 
           (double)(end_time - start_time) / CLOCKS_PER_SEC);
    
    return 0;
}