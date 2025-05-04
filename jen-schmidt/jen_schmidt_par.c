/*
 * To compile this file with OpenMP support, use:
 * gcc -O3 -fopenmp jen_schmidt_par.c -o jen_schmidt_par && ./jen_schmidt_par
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

// Stack implementation for DFS
typedef struct {
    int* array;
    int top;
    int capacity;
    omp_lock_t lock;  // Lock for thread safety
} Stack;

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
    omp_lock_t lock;  // Lock for thread safety
} ComponentList;

// Initialize a stack with thread safety
void init_stack(Stack* stack, int capacity) {
    stack->array = (int*)malloc(capacity * sizeof(int));
    if (!stack->array) {
        fprintf(stderr, "Memory allocation failed for stack\n");
        exit(EXIT_FAILURE);
    }
    stack->capacity = capacity;
    stack->top = -1;
    omp_init_lock(&stack->lock);
}

// Check if the stack is empty
bool is_stack_empty(Stack* stack) {
    return stack->top == -1;
}

// Push an item to stack with thread safety
void push(Stack* stack, int item) {
    omp_set_lock(&stack->lock);
    if (stack->top == stack->capacity - 1) {
        stack->capacity *= 2;
        stack->array = (int*)realloc(stack->array, stack->capacity * sizeof(int));
        if (!stack->array) {
            fprintf(stderr, "Memory allocation failed for stack resize\n");
            exit(EXIT_FAILURE);
        }
    }
    stack->array[++stack->top] = item;
    omp_unset_lock(&stack->lock);
}

// Pop an item from stack with thread safety
int pop(Stack* stack) {
    omp_set_lock(&stack->lock);
    if (is_stack_empty(stack)) {
        fprintf(stderr, "Stack underflow\n");
        exit(EXIT_FAILURE);
    }
    int result = stack->array[stack->top--];
    omp_unset_lock(&stack->lock);
    return result;
}

// Free stack memory
void free_stack(Stack* stack) {
    free(stack->array);
    stack->array = NULL;
    stack->top = -1;
    stack->capacity = 0;
    omp_destroy_lock(&stack->lock);
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

// Initialize component list with thread safety
void init_component_list(ComponentList* list, int initial_capacity) {
    list->count = 0;
    list->capacity = initial_capacity;
    list->components = (BiconnectedComponent*)malloc(initial_capacity * sizeof(BiconnectedComponent));
    if (!list->components) {
        fprintf(stderr, "Failed to allocate component list\n");
        exit(EXIT_FAILURE);
    }
    omp_init_lock(&list->lock);
}

// Add a component to the component list with thread safety
void add_component_to_list(ComponentList* list, BiconnectedComponent* component) {
    omp_set_lock(&list->lock);
    if (list->count == list->capacity) {
        list->capacity *= 2;
        list->components = (BiconnectedComponent*)realloc(list->components, 
                                                         list->capacity * sizeof(BiconnectedComponent));
        if (!list->components) {
            fprintf(stderr, "Failed to reallocate component list\n");
            exit(EXIT_FAILURE);
        }
    }
    
    init_component(&list->components[list->count], component->size);
    for (int i = 0; i < component->size; i++) {
        add_edge_to_component(&list->components[list->count], component->edges[i]);
    }
    list->count++;
    omp_unset_lock(&list->lock);
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
    omp_destroy_lock(&list->lock);
}

// Thread-local function for finding biconnected components
void jen_schmidt_dfs(int u, int* disc, int* low, int* parent, bool* visited, 
                    bool* is_articulation, Edge* edges, Stack* edge_stack, 
                    int* time_ptr, Graph* graph, ComponentList* result) {
    int children = 0;
    visited[u] = true;
    
    // Atomically get and increment time
    int discovery_time;
    #pragma omp atomic capture
        discovery_time = (*time_ptr)++;
    
    disc[u] = low[u] = discovery_time;
    
    // Process all adjacent vertices
    AdjListNode* adj = graph->adjacency_list[u];
    
    while (adj) {
        int v = adj->vertex;
        int edge_id = adj->edge_id;
        
        // Store edge information if not already processed
        if (edges[edge_id].src == 0 && edges[edge_id].dest == 0) {
            edges[edge_id].src = u;
            edges[edge_id].dest = v;
        }
        
        // If v is not visited yet, make it a child of u in DFS tree
        if (!visited[v]) {
            children++;
            parent[v] = u;
            
            // Store the edge in the stack
            push(edge_stack, edge_id);
            
            // Recursive DFS
            jen_schmidt_dfs(v, disc, low, parent, visited, is_articulation, 
                          edges, edge_stack, time_ptr, graph, result);
            
            // Check if subtree rooted with v has a connection to
            // one of the ancestors of u
            if (low[v] >= disc[u]) {
                is_articulation[u] = true;
                
                // Create a new biconnected component
                BiconnectedComponent component;
                init_component(&component, 10);
                
                // Add current edge to component
                add_edge_to_component(&component, edge_id);
                
                // Pop edges from stack until we get back to the current edge
                int e;
                do {
                    omp_set_lock(&edge_stack->lock);
                    e = edge_stack->array[edge_stack->top--];
                    omp_unset_lock(&edge_stack->lock);
                    
                    add_edge_to_component(&component, e);
                } while (e != edge_id);
                
                // Add component to the result
                add_component_to_list(result, &component);
                free(component.edges);
            }
            
            // Update low value of u
            if (low[v] < low[u]) {
                low[u] = low[v];
            }
        } 
        // If v is already visited and not the parent of u
        else if (v != parent[u]) {
            // Update low value of u
            if (disc[v] < low[u]) {
                low[u] = disc[v];
            }
            
            // If this is a back edge, add it to the stack
            if (disc[v] < disc[u]) {
                push(edge_stack, edge_id);
            }
        }
        
        adj = adj->next;
    }
    
    // Root of DFS tree is an articulation point if it has more than one child
    if (parent[u] == -1 && children > 1) {
        is_articulation[u] = true;
    }
}

// Parallel version of Jen-Schmidt algorithm
void jen_schmidt_biconnected_components_parallel(Graph* graph, ComponentList* result) {
    int n = graph->num_nodes;
    int m = graph->num_edges;
    
    // Shared arrays for all threads
    int* disc = (int*)calloc(n, sizeof(int));
    int* low = (int*)calloc(n, sizeof(int));
    int* parent = (int*)malloc(n * sizeof(int));
    bool* visited = (bool*)calloc(n, sizeof(bool));
    bool* is_articulation = (bool*)calloc(n, sizeof(bool));
    Edge* edges = (Edge*)calloc(m, sizeof(Edge));
    
    if (!disc || !low || !parent || !visited || !is_articulation || !edges) {
        fprintf(stderr, "Memory allocation failed for arrays\n");
        exit(EXIT_FAILURE);
    }
    
    // Initialize parent array
    #pragma omp parallel for
    for (int i = 0; i < n; i++) {
        parent[i] = -1;
    }
    
    // Global time counter
    int time = 0;
    
    // Identify connected components
    int* component_ids = (int*)malloc(n * sizeof(int));
    int component_count = 0;
    
    if (!component_ids) {
        fprintf(stderr, "Memory allocation failed for component_ids\n");
        exit(EXIT_FAILURE);
    }
    
    // Initialize all nodes to unassigned component
    #pragma omp parallel for
    for (int i = 0; i < n; i++) {
        component_ids[i] = -1;
    }
    
    // Find connected components using BFS
    for (int i = 0; i < n; i++) {
        if (component_ids[i] == -1) {
            int component_id = component_count++;
            
            // BFS to identify all nodes in this component
            int* queue = (int*)malloc(n * sizeof(int));
            int front = 0, back = 0;
            
            if (!queue) {
                fprintf(stderr, "Memory allocation failed for queue\n");
                exit(EXIT_FAILURE);
            }
            
            queue[back++] = i;
            component_ids[i] = component_id;
            
            while (front != back) {
                int u = queue[front++];
                
                AdjListNode* temp = graph->adjacency_list[u];
                while (temp) {
                    int v = temp->vertex;
                    if (component_ids[v] == -1) {
                        queue[back++] = v;
                        component_ids[v] = component_id;
                    }
                    temp = temp->next;
                }
            }
            
            free(queue);
        }
    }
    
    printf("Found %d connected components\n", component_count);
    
    // Process each connected component in parallel
    #pragma omp parallel
    {
        #pragma omp single
        {
            for (int c = 0; c < component_count; c++) {
                #pragma omp task
                {
                    // Find a starting node for this component
                    int start_node = -1;
                    for (int i = 0; i < n; i++) {
                        if (component_ids[i] == c) {
                            start_node = i;
                            break;
                        }
                    }
                    
                    if (start_node != -1) {
                        // Create thread-local stack
                        Stack edge_stack;
                        init_stack(&edge_stack, m);
                        
                        // Run DFS for this connected component
                        jen_schmidt_dfs(start_node, disc, low, parent, visited, is_articulation,
                                      edges, &edge_stack, &time, graph, result);
                        
                        free_stack(&edge_stack);
                    }
                }
            }
        }
    }
    
    // Clean up
    free(disc);
    free(low);
    free(parent);
    free(visited);
    free(is_articulation);
    free(component_ids);
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
    
    // Print each component in a standardized format matching slota_madduri
    for (int i = 0; i < components->count; i++) {
        printf("Component %d: Size: %d, Nodes: ", i + 1, components->components[i].size);
        
        // Create an array to track which nodes are in this component
        bool* node_in_component = (bool*)calloc(1000, sizeof(bool)); // Assuming max node id < 1000
        if (!node_in_component) {
            fprintf(stderr, "Memory allocation failed for node tracking\n");
            return;
        }
        
        // Mark all nodes that appear in this component
        for (int j = 0; j < components->components[i].size; j++) {
            int edge_id = components->components[i].edges[j];
            int src = edges[edge_id].src;
            int dest = edges[edge_id].dest;
            node_in_component[src] = true;
            node_in_component[dest] = true;
        }
        
        // Print nodes in order
        for (int node = 0; node < 1000; node++) {
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
    
    // Prepare edge array for the results
    Edge* edges = (Edge*)calloc(graph->num_edges, sizeof(Edge));
    if (!edges) {
        fprintf(stderr, "Memory allocation failed for edges array\n");
        free_graph(graph);
        return 1;
    }
    
    ComponentList components;
    printf("Running parallel Jen-Schmidt computation with %d threads...\n", omp_get_max_threads());
    
    init_component_list(&components, 1000);
    
    clock_t compute_start = clock();
    jen_schmidt_biconnected_components_parallel(graph, &components);
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