import networkx as nx
import matplotlib.pyplot as plt
import itertools
import collections
import random
import pickle

MOTIF = [
    nx.Graph([(1,2),(2,3),(3,1),(1,0)]),
    #nx.Graph([(1,2),(1,3),(1,0)]),
    nx.Graph([(1,2),(1,0),(3,0)]),
    nx.complete_graph(4),
    nx.Graph([(1,2),(2,3),(3,0),(1,0),(2,4)])
    #nx.Graph([(1,2),(2,3),(3,0),(1,0)])
]



def run(path):
    result = collections.defaultdict(set)
    result = []


    g = nx.read_edgelist(path)
    nodes = g.nodes()
    for c, n in enumerate(itertools.combinations(nodes, 4)):
        if c % 1000 == 0:
            print path, c
        subg = g.subgraph(n)
        if not nx.is_connected(subg):
            continue
        for index, motif in enumerate(MOTIF):
            if nx.is_isomorphic(subg, motif):
                result.append(subg.nodes())

    print len(result), "and done."
    f = open(path + ".pkl", "wb")
    pickle.dump(result, f)
    f.close()

#run(path = "./BA1.txt")
#run(path = "./BA2.txt")

def get_graph_and_motif(path):
    g = nx.read_edgelist(path)
    motif = pickle.load(open(path + ".pkl", "rb"))
    print "graph", path, len(motif), "motifs loaded."
    return g, motif

def get_node_motif_degree(motif):
    motif_score = collections.defaultdict(int)
    for m in motif:
        for n in m:
            motif_score[int(n)] += 1
    return motif_score

def score(path):
    g, motif = get_graph_and_motif(path)
    return get_node_motif_degree(motif)

s1 = score(path = "./BA1.txt")
s2 = score(path = "./BA2.txt")
ps1 = plt.subplot(211)
ps1.plot(s1.keys(), s1.values())
ps2 = plt.subplot(212)
ps2.plot(s2.keys(), s2.values())
plt.show()

def random_revealed(g, alpha=0.8):
    nodes = g.nodes()
    return set(random.sample(nodes, int(alpha*len(nodes))))

def flood_by_motif(g, motif, motif_score, flood_score):

    for m in motif:
        revealed = [n for n in m if n not in unmatch_nodes]
        unmatched = set(m) - set(revealed)
        s = sum([1.0/motif_score[n] for n in revealed])
        for un in unmatch_nodes:
            flood_score[un] += s



def predict_unmatch(path1, path2):
    g1, motif1 = get_graph_and_motif(path1)
    motif_score1 = get_node_motif_degree(motif1)
    g2, motif2 = get_graph_and_motif(path2)
    motif_score2 = get_node_motif_degree(motif2)

    revealed_nodes = random_revealed(g1)
    unmatch_nodes = set(g1.nodes()) - revealed_nodes

    for n in unmatch_nodes:

