/**
 * Implementation of the Slota-Madduri algorithm for finding articulation points and 
 * biconnected components in graphs.
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

// Structure to represent a biconnected component
typedef struct {
    int* nodes;
    int size;
    int capacity;
} BiconnectedComponent;

// Structure to hold multiple biconnected components
typedef struct {
    BiconnectedComponent* components;
    int count;
    int capacity;
} BiconnectedComponentList;

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

// Initialize a biconnected component
void init_biconnected_component(BiconnectedComponent *component, int initial_capacity) {
    component->size = 0;
    component->capacity = initial_capacity;
    component->nodes = (int*)malloc(initial_capacity * sizeof(int));
    if (!component->nodes) {
        fprintf(stderr, "Failed to allocate component nodes\n");
        exit(EXIT_FAILURE);
    }
}

// Add a node to a biconnected component
void add_node_to_component(BiconnectedComponent *component, int node) {
    if (component->size == component->capacity) {
        component->capacity *= 2;
        component->nodes = (int*)realloc(component->nodes, component->capacity * sizeof(int));
        if (!component->nodes) {
            fprintf(stderr, "Failed to reallocate component nodes\n");
            exit(EXIT_FAILURE);
        }
    }
    
    // Check if node is already in component
    for (int i = 0; i < component->size; i++) {
        if (component->nodes[i] == node) return;
    }
    
    component->nodes[component->size++] = node;
}

// Free component memory
void free_biconnected_component(BiconnectedComponent *component) {
    free(component->nodes);
    component->nodes = NULL;
    component->size = 0;
    component->capacity = 0;
}

// Initialize component list
void init_component_list(BiconnectedComponentList *list, int initial_capacity) {
    list->count = 0;
    list->capacity = initial_capacity;
    list->components = (BiconnectedComponent*)malloc(initial_capacity * sizeof(BiconnectedComponent));
    if (!list->components) {
        fprintf(stderr, "Failed to allocate component list\n");
        exit(EXIT_FAILURE);
    }
}

// Add a component to the component list
void add_component_to_list(BiconnectedComponentList *list, BiconnectedComponent *component) {
    if (list->count == list->capacity) {
        list->capacity *= 2;
        list->components = (BiconnectedComponent*)realloc(list->components, list->capacity * sizeof(BiconnectedComponent));
        if (!list->components) {
            fprintf(stderr, "Failed to reallocate component list\n");
            exit(EXIT_FAILURE);
        }
    }
    
    // Create a copy of the component
    init_biconnected_component(&list->components[list->count], component->size);
    for (int i = 0; i < component->size; i++) {
        add_node_to_component(&list->components[list->count], component->nodes[i]);
    }
    list->count++;
}

// Free component list memory
void free_component_list(BiconnectedComponentList *list) {
    for (int i = 0; i < list->count; i++) {
        free_biconnected_component(&list->components[i]);
    }
    free(list->components);
    list->components = NULL;
    list->count = 0;
    list->capacity = 0;
}

// Slota-Madduri algorithm to find articulation points and biconnected components
void find_articulation_points(Graph *graph, bool *is_articulation_point, BiconnectedComponentList *components) {
    int *discovery_time = (int*)malloc(graph->num_nodes * sizeof(int));
    int *low_value = (int*)malloc(graph->num_nodes * sizeof(int));
    int *parent = (int*)malloc(graph->num_nodes * sizeof(int));
    bool *visited = (bool*)calloc(graph->num_nodes, sizeof(bool));
    
    if (!discovery_time || !low_value || !parent || !visited) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }
    
    // Initialize arrays
    for (int i = 0; i < graph->num_nodes; i++) {
        parent[i] = -1

