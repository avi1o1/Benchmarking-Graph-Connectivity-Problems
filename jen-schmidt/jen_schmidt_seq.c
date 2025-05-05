/**
 * Implementation of the Jen-Schmidt algorithm for finding biconnected components in graphs.
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

// Stack implementation for DFS
typedef struct {
    int* array;
    int top;
    int capacity;
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
} ComponentList;

// Initialize a stack
void init_stack(Stack* stack, int capacity) {
    stack->array = (int*)malloc(capacity * sizeof(int));
    if (!stack->array) {
        fprintf(stderr, "Memory allocation failed for stack\n");
        exit(EXIT_FAILURE);
    }
    stack->capacity = capacity;
    stack->top = -1;
}

// Check if the stack is empty
bool is_stack_empty(Stack* stack) {
    return stack->top == -1;
}

// Push an item to stack
void push(Stack* stack, int item) {
    if (stack->top == stack->capacity - 1) {
        stack->capacity *= 2;
        stack->array = (int*)realloc(stack->array, stack->capacity * sizeof(int));
        if (!stack->array) {
            fprintf(stderr, "Memory allocation failed for stack resize\n");
            exit(EXIT_FAILURE);
        }
    }
    stack->array[++stack->top] = item;
}

// Pop an item from stack
int pop(Stack* stack) {
    if (is_stack_empty(stack)) {
        fprintf(stderr, "Stack underflow\n");
        exit(EXIT_FAILURE);
    }
    return stack->array[stack->top--];
}

// Free stack memory
void free_stack(Stack* stack) {
    free(stack->array);
    stack->array = NULL;
    stack->top = -1;
    stack->capacity = 0;
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
void add_component_to_list(ComponentList* list, BiconnectedComponent* component) {
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

// Jen-Schmidt algorithm for finding biconnected components
void jen_schmidt_biconnected_components(Graph* graph, ComponentList* result) {
    int n = graph->num_nodes;
    int m = graph->num_edges;
    
    // Arrays for DFS traversal
    int* disc = (int*)malloc(n * sizeof(int));        // Discovery time
    int* low = (int*)malloc(n * sizeof(int));         // Earliest reachable ancestor
    int* parent = (int*)malloc(n * sizeof(int));      // Parent in DFS tree
    bool* visited = (bool*)calloc(n, sizeof(bool));   // Visited nodes
    bool* is_articulation = (bool*)calloc(n, sizeof(bool)); // Articulation points
    
    // Edge information for result
    Edge* edges = (Edge*)malloc(m * sizeof(Edge));
    
    // Stack for tracking edges in current biconnected component
    Stack edge_stack;
    init_stack(&edge_stack, m);
    
    if (!disc || !low || !parent || !visited || !is_articulation || !edges) {
        fprintf(stderr, "Memory allocation failed for DFS arrays\n");
        exit(EXIT_FAILURE);
    }
    
    for (int i = 0; i < n; i++) {
        parent[i] = -1;
    }
    
    // Initialize result list
    init_component_list(result, 10);  // Start with reasonable capacity
    
    int time = 0;
    
    // Helper function to process biconnected components (implemented as nested function)
    void process_component(int u, int v, int edge_id) {
        BiconnectedComponent component;
        init_component(&component, 10);
        
        // Add this edge to the component
        add_edge_to_component(&component, edge_id);
        
        // Pop edges from stack until we get back to the current edge
        int e;
        do {
            e = pop(&edge_stack);
            add_edge_to_component(&component, e);
        } while (e != edge_id);
        
        // Add component to the result
        add_component_to_list(result, &component);
        free(component.edges);
    }
    
    // DFS function for finding biconnected components (implemented as nested function)
    void biconnected_dfs(int u) {
        visited[u] = true;
        disc[u] = low[u] = ++time;
        
        int children = 0;  // Count of children in DFS tree
        
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
                push(&edge_stack, edge_id);
                
                // Recursive DFS
                biconnected_dfs(v);
                
                // Check if subtree rooted with v has a connection to
                // one of the ancestors of u
                if (low[v] >= disc[u]) {
                    is_articulation[u] = true;
                    process_component(u, v, edge_id);
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
                    push(&edge_stack, edge_id);
                }
            }
            
            adj = adj->next;
        }
        
        // Root of DFS tree is an articulation point if it has more than one child
        if (parent[u] == -1 && children > 1) {
            is_articulation[u] = true;
        }
    }
    
    // Process each connected component
    for (int i = 0; i < n; i++) {
        if (!visited[i]) {
            biconnected_dfs(i);
        }
    }
    
    // Clean up
    free(disc);
    free(low);
    free(parent);
    free(visited);
    free(is_articulation);
    free_stack(&edge_stack);
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
    
    // Initialize edges
    for (int i = 0; i < graph->num_edges; i++) {
        edges[i].src = edges[i].dest = 0;
    }
    
    ComponentList components;
    printf("Running sequential Jen-Schmidt computation...\n");
    
    clock_t compute_start = clock();
    jen_schmidt_biconnected_components(graph, &components);
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