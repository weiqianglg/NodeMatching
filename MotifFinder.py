import networkx as nx
import matplotlib.pyplot as plt
import itertools
import collections
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

def score(path):
    g = nx.read_edgelist(path)
    motif = pickle.load(open(path + ".pkl", "rb"))
    print len(motif), "motifs loaded."
    motif_score = collections.defaultdict(int)
    for m in motif:
        for n in m:
            motif_score[int(n)] += 1
    return motif_score

s1 = score(path = "./BA1.txt")
s2 = score(path = "./BA2.txt")
ps1 = plt.subplot(211)
ps1.plot(s1.keys(), s1.values())
ps2 = plt.subplot(212)
ps2.plot(s2.keys(), s2.values())
plt.show()

