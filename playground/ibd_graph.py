import os
import pandas as pd
import networkx as nx
from infomap import Infomap
import matplotlib.pyplot as plt
import seaborn as sns

sns.set(style="whitegrid")
lpa_dir=os.environ.get("LPA")

match_file = f"{lpa_dir}/ilash/ilash_ch6/output.match"

# Define the correct columns for iLASH output
cols = [
    "id1", "hap1", "id2", "hap2",
    "chr", "start", "end",
    "genomic_start", "genomic_end",
    "cm", "snp_count"
]

# Read with proper whitespace delimiter
df = pd.read_csv(match_file, delim_whitespace=True, names=cols)

# Initialize an undirected graph
G = nx.Graph()

# Add edges with weights (e.g., total shared cM)
for _, row in df.iterrows():
    u, v, cm = row["id1"], row["id2"], row["cm"]
    if G.has_edge(u, v):
        G[u][v]["weight"] += cm
    else:
        G.add_edge(u, v, weight=cm)


# Step 1: Relabel nodes with integers for Infomap
node_to_int = {node: idx for idx, node in enumerate(G.nodes())}
int_to_node = {idx: node for node, idx in node_to_int.items()}

# Step 2: Create new graph with integer node IDs
G_int = nx.relabel_nodes(G, node_to_int)

# Step 3: Run Infomap on the integer-labeled graph
im = Infomap()

for u, v, data in G_int.edges(data=True):
    im.add_link(u, v, float(data.get("weight", 1.0)))

im.run()

# Step 4: Extract community assignments and map back to original labels
communities = {int_to_node[node.node_id]: node.module_id for node in im.nodes}

# Step 5: Assign community attribute to original graph
nx.set_node_attributes(G, communities, "community")

# Optional: Check how many communities were detected
print(f"Number of communities detected: {len(set(communities.values()))}")

# Get communities and layout
colors = [G.nodes[n].get("community", 0) for n in G.nodes()]

pos = nx.spring_layout(G, seed=42)
# pos = nx.kamada_kawai_layout(G, weight='weight')
# pos = nx.forceatlas2_layout(G, weight='weight')
# pos = nx.arf_layout(G, weight='weight')
# pos = nx.spectral_layout(G, weight='weight')
# pos = nx.circular_layout(G)
# pos = nx.random_layout(G, seed=42)

plt.figure(figsize=(12, 12))
nx.draw_networkx_nodes(G, pos, node_size=2, node_color=colors, cmap="tab20")
plt.title("IBD Network with Communities (InfoMap)", fontsize=16)
plt.axis("off")
plt.savefig(f"{lpa_dir}/playground/figures/ibd_network_infomap_spring.png", dpi=300, bbox_inches='tight')
