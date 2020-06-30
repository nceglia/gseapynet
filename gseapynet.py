import pandas
import logging
import gseapy as gp
import networkx as nx
import collections
from itertools import combinations
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.gridspec as gridspec
import numpy
import sys
from gseapy.plot import barplot
import argparse

### Change color map
cmap = matplotlib.cm.coolwarm



parser = argparse.ArgumentParser(description='Run ranked GSEA with network output.')
parser.add_argument('--csv', type=str, help='CSV with DE genes.')
parser.add_argument('--sig', type=float, help='Threshold for signficance in DE input.')
parser.add_argument('--sig-col', type=str, help='Name of column with signficance in input csv.')
parser.add_argument('--minfc', type=float, help='Threshold for fold change column in csv input.')
parser.add_argument('--minfc-col', type=str, help='Name of column with fold change in input csv.')
parser.add_argument('--pathway-fdr', type=float, help='FDR threshold for GSEA pathway result.')
parser.add_argument('--png', type=str, help='Name of png output file. (end in .svg for svg result).')
parser.add_argument('--gmt', type=str, help='Either a custom gmt file or the name of an enrichr library (https://amp.pharm.mssm.edu/Enrichr/#stats)')
args = parser.parse_args()


plt.rcParams['svg.fonttype'] = 'none'

logger = logging.getLogger("Pathway Analysis")
logger.setLevel(logging.INFO)
logging.basicConfig(level=logging.INFO)

csv = args.csv
fdr = args.pathway_fdr
minfc = args.minfc
adjp = args.minfc_col
minedge = 1
png = args.png

logger.info("Reading CSV file")
deg = pandas.read_csv(csv)

deg = deg[abs(deg[args.sig_col]) < float(args.sig)]
fig = plt.figure(figsize=(10,7))
gs = gridspec.GridSpec(nrows=2, ncols=2, width_ratios=[30,1], height_ratios=[30, 1])
axes = []
scales = []
axes.append(fig.add_subplot(gs[0, 0]))
axes.append(fig.add_subplot(gs[1, 0]))
axes.append(fig.add_subplot(gs[0, 1]))
margin=0.05
fig.subplots_adjust(margin, margin, 1.-margin, 1.-margin)

deg = deg[abs(deg[args.minfc_col]) > float(minfc)]
data = deg[["gene",args.minfc_col]]

data = data.sort_values(args.minfc_col)

pre_res = gp.prerank(rnk=data, gene_sets=args.gmt,
                    processes=4,
                    permutation_num=100,
                    outdir="prerank", format='png', seed=6,min_size=1,max_size=30000000000)


print(pre_res.res2d.sort_index().head())
result = pre_res.res2d
result = result[result["fdr"]<float(args.pathway_fdr)] #### ADD ME BACK
#
result.to_csv(png.replace(".png",".csv").replace(".svg",".csv"),header=True)
pathways = result.index
nodes = collections.defaultdict(dict)
edges = dict()

minnes = min(result["nes"])
maxnes = max(result["nes"])

for term in pathways:
    hallmark = term
    print(result.loc[hallmark]["genes"])
    term = term.replace("HALLMARK_","").replace("_"," ")
    nodes[term]["genes"] = result.loc[hallmark]["genes"].split(";")
    nodes[term]["NES"] = result.loc[hallmark]["nes"]
G = nx.Graph()
node_color = []
node_order = []
node_size = []
for term in nodes:
    node_order.append(term)
    node_color.append(float(nodes[term]["NES"]))
    node_size.append(abs(float(nodes[term]["NES"] * 500.0)))
    G.add_node(term, nes=float(nodes[term]["NES"]))

if len(node_color) == 1:
    print("No Enriched Pathways! Try relaxing the paramters.")
    exit(0)

edge_order = []
edge_color = []
overlaps = []
for edge in combinations(node_order,2):
    overlap = set(nodes[edge[0]]["genes"]).intersection(nodes[edge[1]]["genes"])
    if len(overlap) > 0:
        overlaps.append(len(overlap))
        G.add_edge(edge[0], edge[1], weight=len(overlap))
        print(edge[0],edge[1],len(overlap), overlap)
        edge_color.append(len(overlap))
        edge_order.append(edge)

G.remove_nodes_from(list(nx.isolates(G)))
if len(overlaps) == 0:
    maxvedge = 1
else:
    maxvedge = max(overlaps)
_node_order = []
_node_size = []
_node_color = []
for nord, nsiz, ncol in zip(node_order, node_size, node_color):
    if nord in G.nodes():
        _node_order.append(nord)
        _node_size.append(nsiz)
        _node_color.append(ncol)
node_order = _node_order
node_size = _node_size
node_color = _node_color
print(len(node_order), len(node_size), len(G.nodes()))
pos=nx.spring_layout(G,k=0.1,scale=0.03,fixed=None)
axes[0].set_title(args.gmt, fontsize=13)
axes[0].axis('equal')

nx.draw(G,pos,ax=axes[0], nodelist=node_order, node_size=node_size, vmin=minnes, vmax=maxnes,edgelist=edge_order, node_color=node_color, edge_color=edge_color, cmap=cmap, edge_cmap=plt.cm.Greys, edge_vmin=0, edge_vmax=maxvedge, with_labels=True, width=2,font_size=11)


norm = matplotlib.colors.Normalize(vmin=-8, vmax=18)
cb1 = matplotlib.colorbar.ColorbarBase(axes[1], cmap=cmap,
                                        norm=norm,
                                        orientation='horizontal')
axes[1].set_title("Normalized Enrichment Score",fontsize=13)
if overlaps == []:
    overlaps.append(12)

cmap2 = matplotlib.cm.binary
norm2 = matplotlib.colors.Normalize(vmin=0, vmax=max(overlaps))
cb2 = matplotlib.colorbar.ColorbarBase(axes[2], cmap=cmap2,
                                        norm=norm2,
                                        orientation='vertical', drawedges=False)
cb2.set_ticks(list(range(0,maxvedge,1)))
axes[2].set_title("Common Genes",fontsize=10)
plt.tight_layout()
plt.savefig(png)