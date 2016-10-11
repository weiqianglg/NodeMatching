import os.path
import networkx as nx


def get_store_filename(path):
    dirn, filen = os.path.split(path)
    new_filen = filen[:-len(".txt")] + "-csub.txt"
    return os.path.join(dirn, new_filen)


gpath1 = r"./multiplex_embeddings_data/Drosophila Melanogaster/l_1_giant-links.txt"
gpath2 = r"./multiplex_embeddings_data/Drosophila Melanogaster/l_2_giant-links.txt"
g1 = nx.read_edgelist(gpath1)
g2 = nx.read_edgelist(gpath2)
common_nodes = set(g1.nodes()) & set(g2.nodes())
sub_g1 = g1.subgraph(common_nodes)
sub_g2 = g2.subgraph(common_nodes)
print sub_g1.number_of_nodes(), sub_g2.number_of_nodes()
nx.write_edgelist(sub_g1, get_store_filename(gpath1))
nx.write_edgelist(sub_g2, get_store_filename(gpath2))
print "done."
