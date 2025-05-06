#!/usr/bin/env python3
# filepath: /home/avilol/Downloads/GitHub/Benchmarking-Graph-Connectivity-Problems/tmp.py
# Directed graph (each unordered pair of nodes is saved once): CA-HepPh.txt 
# Collaboration network of Arxiv High Energy Physics category (there is an edge if authors coauthored at least one paper)
# Nodes: 12008 Edges: 237010
# FromNodeId	ToNodeId

import os

def renumber_graph_nodes(input_file, output_file):
    """
    Renumber the nodes in a graph from arbitrary IDs to sequential IDs (0 to N-1)
    
    Args:
        input_file (str): Path to the input graph file
        output_file (str): Path to save the renumbered graph
    """
    # Read the edges from the file
    edges = []
    unique_nodes = set()
    
    print(f"Reading edges from {input_file}...")
    with open(input_file, 'r') as f:
        for line in f:
            # Skip empty lines
            if not line.strip():
                continue
                
            # Parse edge
            try:
                source, target = map(int, line.strip().split())
                edges.append((source, target))
                unique_nodes.add(source)
                unique_nodes.add(target)
            except ValueError:
                # Skip lines that don't have exactly two integers
                continue
    
    # Create a mapping from original node IDs to sequential IDs
    print(f"Creating node mapping for {len(unique_nodes)} unique nodes...")
    node_mapping = {original_id: new_id for new_id, original_id in enumerate(sorted(unique_nodes))}
    
    # Renumber the edges
    print("Renumbering edges...")
    new_edges = [(node_mapping[src], node_mapping[dst]) for src, dst in edges]
    
    # Write the renumbered edges to the output file
    print(f"Writing renumbered graph to {output_file}...")
    with open(output_file, 'w') as f:
        for src, dst in new_edges:
            f.write(f"{src} {dst}\n")
    
    print(f"Done! Renumbered {len(edges)} edges and {len(unique_nodes)} nodes.")
    print(f"New node IDs range from 0 to {len(unique_nodes) - 1}")

if __name__ == "__main__":
    input_file = "./datasets/CA-HepPh.txt"
    output_file = "./datasets/CA-HepPh-renumbered.txt"
    
    # Ensure input file exists
    if not os.path.exists(input_file):
        print(f"Error: Input file {input_file} does not exist.")
        exit(1)
    
    renumber_graph_nodes(input_file, output_file)