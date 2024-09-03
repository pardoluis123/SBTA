import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from Data_Manip import Data_Manipulator
import scipy.sparse as sp
from Convenience import n2_residue_numbers, restrained_residue_list
import mdtraj as md

# Load the arrays
array_CCUCGU = sp.load_npz("/zfshomes/lperez/fingerprint/H_Print/CCUCGU_Replicate_Average.npz")
array_CCUGCU = sp.load_npz("/zfshomes/lperez/fingerprint/H_Print/CCUGCU_Replicate_Average.npz")


# Process the arrays using Data_Manipulator
CCUGCU = Data_Manipulator(
    array_one=array_CCUGCU,
    array_two=array_CCUCGU,
    residues_to_filter=n2_residue_numbers,
    res_filtered=restrained_residue_list,
    topology=md.load_topology("/home66/kscopino/AMBER22/CODONS/CCUGCU_G34/TLEAP/5JUP_N2_GCU_nowat.prmtop")
)

CCUGCU_For_Network = CCUGCU.filter_array()
CCUGCU.array_two=CCUGCU.process_input(array=CCUGCU.array_two)
CCUCGU_For_Network = CCUGCU.filter_array(array=CCUGCU.array_two)
CCU_network=np.copy(CCUGCU_For_Network)

CCU_network[1:,1:]=CCUGCU_For_Network[1:,1:] - CCUCGU_For_Network[1:,1:]

print(CCU_network[0,:])

# Remove the first row and the first column
CCU_network = np.delete(CCUGCU_For_Network, 0, axis=0)  # Remove the first row
CCU_network = np.delete(CCUGCU_For_Network, 0, axis=1)  # Remove the first column

# Initialize an empty edge list
edges = []

# Get the number of nodes
num_nodes = CCU_network.shape[0]

# Iterate over the matrix to extract edges
for i in range(num_nodes-1):
    for j in range(num_nodes-1):
        if i != j:  # Assuming no self-loops
            weight = CCU_network[i, j]
            if weight != 0:
                edges.append((i + 1, j + 1, weight))  # Renumber nodes from 0-based to 1-based

# Initialize an undirected graph
G = nx.Graph()

# Add edges to the graph
for edge in edges:
    G.add_edge(edge[0], edge[1], weight=edge[2])

# Draw the graph with edge weights
pos = nx.spring_layout(G)  # You can use other layouts like circular_layout, shell_layout, etc.

plt.figure(figsize=(8, 6))

# Get the edge weights
weights = [G[u][v]['weight'] for u, v in G.edges()]

# Normalize the weights to use them for edge width
max_weight = max(weights)
min_weight = min(weights)
normalized_weights = [(weight - min_weight) / (max_weight - min_weight) * 5 + 1 for weight in weights]

nx.draw(G, pos, with_labels=True, node_size=700, node_color='lightblue', font_size=10, font_weight='bold',
        width=normalized_weights)
plt.title('Pairwise Comparison Network')
plt.savefig('difference_Network_Hbonds_renumbered_undirected.png', dpi=300, bbox_inches='tight')
plt.clf()
