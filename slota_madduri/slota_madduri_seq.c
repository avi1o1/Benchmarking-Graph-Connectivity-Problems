/**
 * Implementation of the Slota-Madduri algorithm for finding maximal cliques in graphs.
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
    struct AdjListNode* next;
} AdjListNode;

// Graph structure with adjacency list
typedef struct {
    int num_nodes;
    AdjListNode** adjacency_list;
    int* degrees;
} Graph;

// Clique structure
typedef struct {
    int* nodes;
    int size;
    int capacity;
} Clique;

// Structure to hold multiple cliques
typedef struct {
    Clique* cliques;
    int count;
    int capacity;
} CliqueList;

// Initialize an empty graph with n nodes
void init_graph(Graph *graph, int n) {
    graph->num_nodes = n;
    graph->adjacency_list = (AdjListNode**)malloc(n * sizeof(AdjListNode*));
    if (!graph->adjacency_list) {
        fprintf(stderr, "Failed to allocate adjacency list\n");
        exit(EXIT_FAILURE);
    }
    
    graph->degrees = (int*)calloc(n, sizeof(int));
    if (!graph->degrees) {
        fprintf(stderr, "Failed to allocate degrees array\n");
        free(graph->adjacency_list);
        exit(EXIT_FAILURE);
    }
    
    for (int i = 0; i < n; i++) {
        graph->adjacency_list[i] = NULL;
    }
}

// Create a new adjacency list node
AdjListNode* new_adj_list_node(int dest) {
    AdjListNode* new_node = (AdjListNode*)malloc(sizeof(AdjListNode));
    if (!new_node) {
        fprintf(stderr, "Failed to allocate adjacency list node\n");
        exit(EXIT_FAILURE);
    }
    new_node->vertex = dest;
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
        
        // Add edge u -> v
        AdjListNode* new_node = new_adj_list_node(v);
        new_node->next = graph->adjacency_list[u];
        graph->adjacency_list[u] = new_node;
        graph->degrees[u]++;
        
        // Add edge v -> u
        new_node = new_adj_list_node(u);
        new_node->next = graph->adjacency_list[v];
        graph->adjacency_list[v] = new_node;
        graph->degrees[v]++;
    }
}

// Check if two nodes are connected
bool has_edge(Graph *graph, int u, int v) {
    AdjListNode* temp = graph->adjacency_list[u];
    while (temp) {
        if (temp->vertex == v) {
            return true;
        }
        temp = temp->next;
    }
    return false;
}

// Initialize an empty clique
void init_clique(Clique *clique, int initial_capacity) {
    clique->size = 0;
    clique->capacity = initial_capacity;
    clique->nodes = (int*)malloc(initial_capacity * sizeof(int));
    if (!clique->nodes) {
        fprintf(stderr, "Failed to allocate clique nodes\n");
        exit(EXIT_FAILURE);
    }
}

// Add a node to a clique
void add_node_to_clique(Clique *clique, int node) {
    if (clique->size == clique->capacity) {
        clique->capacity *= 2;
        clique->nodes = (int*)realloc(clique->nodes, clique->capacity * sizeof(int));
        if (!clique->nodes) {
            fprintf(stderr, "Failed to reallocate clique nodes\n");
            exit(EXIT_FAILURE);
        }
    }
    clique->nodes[clique->size++] = node;
}

// Check if a node is in a clique
bool is_node_in_clique(Clique *clique, int node) {
    for (int i = 0; i < clique->size; i++) {
        if (clique->nodes[i] == node) {
            return true;
        }
    }
    return false;
}

// Free clique memory
void free_clique(Clique *clique) {
    free(clique->nodes);
    clique->nodes = NULL;
    clique->size = 0;
    clique->capacity = 0;
}

// Initialize clique list
void init_clique_list(CliqueList *list, int initial_capacity) {
    list->count = 0;
    list->capacity = initial_capacity;
    list->cliques = (Clique*)malloc(initial_capacity * sizeof(Clique));
    if (!list->cliques) {
        fprintf(stderr, "Failed to allocate clique list\n");
        exit(EXIT_FAILURE);
    }
}

// Add a clique to the clique list
void add_clique_to_list(CliqueList *list, Clique *clique) {
    if (list->count == list->capacity) {
        list->capacity *= 2;
        list->cliques = (Clique*)realloc(list->cliques, list->capacity * sizeof(Clique));
        if (!list->cliques) {
            fprintf(stderr, "Failed to reallocate clique list\n");
            exit(EXIT_FAILURE);
        }
    }
    
    // Create a copy of the clique
    init_clique(&list->cliques[list->count], clique->size);
    for (int i = 0; i < clique->size; i++) {
        add_node_to_clique(&list->cliques[list->count], clique->nodes[i]);
    }
    list->count++;
}

// Free clique list memory
void free_clique_list(CliqueList *list) {
    for (int i = 0; i < list->count; i++) {
        free_clique(&list->cliques[i]);
    }
    free(list->cliques);
    list->cliques = NULL;
    list->count = 0;
    list->capacity = 0;
}

// Check if a clique is maximal
bool is_maximal_clique(Graph *graph, Clique *clique, bool *visited) {
    for (int node = 0; node < graph->num_nodes; node++) {
        if (!is_node_in_clique(clique, node) && !visited[node]) {
            bool can_add = true;
            for (int i = 0; i < clique->size; i++) {
                if (!has_edge(graph, node, clique->nodes[i])) {
                    can_add = false;
                    break;
                }
            }
            if (can_add) {
                return false;  // Found a node that can be added, so not maximal
            }
        }
    }
    return true;
}

// Find a maximal clique starting from a specific node
void find_maximal_clique_from_node(Graph *graph, int start_node, bool *visited, Clique *clique) {
    init_clique(clique, graph->degrees[start_node] + 1);
    add_node_to_clique(clique, start_node);
    
    // Greedily add nodes to form a clique
    for (int node = 0; node < graph->num_nodes; node++) {
        if (node != start_node && !visited[node] && has_edge(graph, start_node, node)) {
            bool can_add = true;
            for (int i = 0; i < clique->size; i++) {
                if (!has_edge(graph, node, clique->nodes[i])) {
                    can_add = false;
                    break;
                }
            }
            if (can_add) {
                add_node_to_clique(clique, node);
            }
        }
    }
    
    // Ensure clique is maximal by trying to add more nodes
    bool expanded;
    do {
        expanded = false;
        for (int node = 0; node < graph->num_nodes; node++) {
            if (!is_node_in_clique(clique, node) && !visited[node]) {
                bool can_add = true;
                for (int i = 0; i < clique->size; i++) {
                    if (!has_edge(graph, node, clique->nodes[i])) {
                        can_add = false;
                        break;
                    }
                }
                if (can_add) {
                    add_node_to_clique(clique, node);
                    expanded = true;
                    break;
                }
            }
        }
    } while (expanded);
}

// Slota-Madduri algorithm for finding maximal cliques - sequential version
void slota_madduri_maximal_cliques(Graph *graph, CliqueList *result) {
    bool *visited = (bool *)calloc(graph->num_nodes, sizeof(bool));
    if (!visited) {
        fprintf(stderr, "Memory allocation failed for visited array\n");
        exit(EXIT_FAILURE);
    }
    
    // Sort nodes by degree (descending)
    int *nodes = (int *)malloc(graph->num_nodes * sizeof(int));
    if (!nodes) {
        fprintf(stderr, "Memory allocation failed for nodes array\n");
        free(visited);
        exit(EXIT_FAILURE);
    }
    
    for (int i = 0; i < graph->num_nodes; i++) {
        nodes[i] = i;
    }
    
    // Sort nodes by degree (descending)
    for (int i = 0; i < graph->num_nodes; i++) {
        for (int j = i + 1; j < graph->num_nodes; j++) {
            if (graph->degrees[nodes[i]] < graph->degrees[nodes[j]]) {
                int temp = nodes[i];
                nodes[i] = nodes[j];
                nodes[j] = temp;
            }
        }
    }
    
    init_clique_list(result, 1000); // Start with reasonable capacity
    
    // Process all nodes
    for (int idx = 0; idx < graph->num_nodes; idx++) {
        int node = nodes[idx];
        
        if (!visited[node]) {
            Clique clique;
            find_maximal_clique_from_node(graph, node, visited, &clique);
            
            if (is_maximal_clique(graph, &clique, visited)) {
                add_clique_to_list(result, &clique);
            }
            
            for (int i = 0; i < clique.size; i++) {
                visited[clique.nodes[i]] = true;
            }
            
            free_clique(&clique);
        }
    }
    
    free(visited);
    free(nodes);
    
    // Post-processing: remove non-maximal cliques
    CliqueList filtered;
    init_clique_list(&filtered, result->count);
    
    for (int i = 0; i < result->count; i++) {
        bool is_subset = false;
        for (int j = 0; j < filtered.count; j++) {
            if (result->cliques[i].size <= filtered.cliques[j].size) {
                int matches = 0;
                for (int k = 0; k < result->cliques[i].size; k++) {
                    for (int l = 0; l < filtered.cliques[j].size; l++) {
                        if (result->cliques[i].nodes[k] == filtered.cliques[j].nodes[l]) {
                            matches++;
                            break;
                        }
                    }
                }
                if (matches == result->cliques[i].size) {
                    is_subset = true;
                    break;
                }
            }
        }
        if (!is_subset) {
            add_clique_to_list(&filtered, &result->cliques[i]);
        }
    }
    
    free_clique_list(result);
    *result = filtered;
}

// Print the cliques in a human-readable format
void print_cliques(CliqueList *cliques) {
    printf("Found %d maximal cliques:\n", cliques->count);
    
    // Sort cliques by size (largest first) before printing
    for (int i = 0; i < cliques->count - 1; i++) {
        for (int j = i + 1; j < cliques->count; j++) {
            if (cliques->cliques[i].size < cliques->cliques[j].size) {
                Clique temp = cliques->cliques[i];
                cliques->cliques[i] = cliques->cliques[j];
                cliques->cliques[j] = temp;
            }
        }
    }
    
    // Print each clique (limit output for very large results)
    int max_to_print = (cliques->count > 100) ? 100 : cliques->count;
    for (int i = 0; i < max_to_print; i++) {
        printf("Clique %d: Size: %d, Nodes: ", i + 1, cliques->cliques[i].size);
        // Print only a few nodes for large cliques
        int nodes_to_print = (cliques->cliques[i].size > 20) ? 20 : cliques->cliques[i].size;
        for (int j = 0; j < nodes_to_print; j++) {
            printf("%d ", cliques->cliques[i].nodes[j]);
        }
        if (nodes_to_print < cliques->cliques[i].size) {
            printf("... (and %d more)", cliques->cliques[i].size - nodes_to_print);
        }
        printf("\n");
    }
    
    if (max_to_print < cliques->count) {
        printf("... (showing only %d out of %d cliques)\n", max_to_print, cliques->count);
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
    free(graph->degrees);
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
    if (fgets(line, sizeof(line), file) == NULL) {  // Skip the first line (number of nodes)
        fprintf(stderr, "Error reading first line\n");
        free(graph);
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
        filename = "datasets/large.txt";
    }
    
    clock_t start_time = clock();
    
    printf("Loading graph from %s...\n", filename);
    Graph* graph = load_graph_from_file(filename);
    
    if (!graph) {
        fprintf(stderr, "Failed to load graph from %s\n", filename);
        return 1;
    }
    
    printf("Graph loaded from %s with %d nodes\n", filename, graph->num_nodes);
    
    CliqueList cliques;
    printf("Running sequential computation...\n");
    
    clock_t compute_start = clock();
    slota_madduri_maximal_cliques(graph, &cliques);
    clock_t compute_end = clock();
    
    printf("Computation time: %.5f seconds\n", 
           (double)(compute_end - compute_start) / CLOCKS_PER_SEC);
    
    print_cliques(&cliques);
    
    // Clean up
    free_clique_list(&cliques);
    free_graph(graph);
    
    clock_t end_time = clock();
    printf("Total execution time: %.5f seconds\n", 
           (double)(end_time - start_time) / CLOCKS_PER_SEC);
    
    return 0;
}

