# -*- coding: utf-8 -*-
import networkx as nx
import random
import itertools
import logging


LAYER_PREFIX = ['0_', '1_', '2_', '3_', '4_', '5_', '']


class MergeMultiGraph(object):

    def __init__(self, *layer_path):
        self.graphs = []
        self.nodes = []
        for index, p in enumerate(layer_path):
            g, g_nodes = self.construct_graph_from_file(p, index)
            self.graphs.append(g)
            self.nodes.append(g_nodes)

    def construct_graph_from_file(self, path, index=-1):
        """read edgelist and rename each node by index."""
        logging.debug("load graph from {}.".format(path))
        g = nx.Graph()
        nodes = set()
        with open(path, "r") as f:
            for line in f.readlines():
                x, y = line.split()
                nodes.add(x)
                nodes.add(y)
                x, y = LAYER_PREFIX[index] + x, LAYER_PREFIX[index] + y
                g.add_edge(x, y)
        logging.info("graph {} has {} nodes and {} edges.".format
                    (index, g.number_of_nodes(), g.number_of_edges()))
        return (g, nodes)

    def random_select_revealed_node(self, alpha, index1, index2):
        """select some nodes as revealed nodes randomly,
        alpha is percent of all same nodes."""
        same_node = self.nodes[index1] & self.nodes[index2]
        s = int(alpha * len(same_node))
        logging.info("graph {}-{} random revealed nodes {}/{}.".format
                     (index1, index2, s, len(same_node)))
        return random.sample(same_node, s)

    def construct_merge_sequence(self):
        """construct merge seq which is a sorted revealed nodes on
        different layer, node in the seq will be merged together."""
        sequence = []
        alpha = 0.5
        g = nx.Graph()
        for i, j in itertools.combinations(list(range(0, len(self.graphs))), 2):
            nodes = self.random_select_revealed_node(alpha, i, j)
            for n in nodes:
                g.add_edge(LAYER_PREFIX[i] + n, LAYER_PREFIX[j] + n)
        while g.size() > 0:
            n = g.nodes()[0]
            all_n = nx.bfs_tree(g, n).nodes()
            sequence.append(
                sorted(all_n, cmp=lambda x, y: cmp(x, y))
                )
            g.remove_nodes_from(all_n)
        return sequence

    def merge_graphs(self):
        """merge revealed nodes to get a new graph, and write to file."""
        g = nx.union_all(self.graphs)
        merge_seq = self.construct_merge_sequence()
        for seq in merge_seq:
            keep, merge = seq[:1], seq[1:]
            for m in merge:
                neighbor_m = g.neighbors(m)
                g.remove_node(m)
                g.add_star(keep + neighbor_m)
            logging.debug("merge nodes {}.".format(seq))
        logging.info("merged graph has {} nodes and {} edges.".format
                    (g.number_of_nodes(), g.number_of_edges()))
        nx.write_edgelist(g, "merged_net.txt", data=False)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    path1 = r"./multiplex_embeddings_data/IPv4_IPv6 Internet/l_1.txt"
    path2 = r"./multiplex_embeddings_data/IPv4_IPv6 Internet/l_2.txt"
    c = MergeMultiGraph(path1, path2)
    c.merge_graphs()






