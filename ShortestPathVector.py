import networkx as nx
import numpy as np
from scipy.spatial.distance import cdist
import logging
import random
from munkres import Munkres


class ShortestPathLengthVector(object):
    def __init__(self, g1, g2):
        self.g1 = nx.read_edgelist(g1)
        self.g2 = nx.read_edgelist(g2)
        logging.info("graph loaded. G1({},{}), G2({},{})".format
                     (self.g1.order(), self.g1.size(), self.g2.order(), self.g2.size()))
        self.same_nodes = set(self.g1.nodes()) & set(self.g2.nodes())
        self.revealed_nodes = self.random_revealed_nodes(0.9)
        self.unmatched_nodes = self.same_nodes - self.revealed_nodes

    def random_revealed_nodes(self, alpha):
        num = int(alpha * len(self.same_nodes))
        logging.info("random revealed nodes {}/{}.".format
                     (num, len(self.same_nodes)))
        return set(random.sample(self.same_nodes, num))

    def cal_spl_array(self, g):
        r = np.zeros((len(self.unmatched_nodes), len(self.revealed_nodes)), dtype=np.uint8)
        for i, n in enumerate(self.unmatched_nodes):
            n2other = nx.single_source_shortest_path_length(g, n)
            for j, revealed in enumerate(self.revealed_nodes):
                if revealed in n2other:
                    r[i, j] = n2other[revealed]
                else:  # no path
                    r[i, j] = 0
        return r

    def cal_correlation(self):
        r1 = self.cal_spl_array(self.g1)
        r2 = self.cal_spl_array(self.g2)
        y = cdist(r1, r2, metric="correlation")
        return y

    def km_match(self):
        distance_matrix = self.cal_correlation()
        total_score, right_score = 0, 0
        right_match = 0
        logging.info("computing min match by KM...")
        km_algo = Munkres()
        indexes = km_algo.compute(distance_matrix.copy())
        for row, column in indexes:
            logging.info("{}-{} with score {}".format(row, column, distance_matrix[row, column]))
            total_score += distance_matrix[row, column]
            right_score += distance_matrix[row, row]
            right_match += 1 if row == column else 0
        logging.info("km score {}, vs right score {}".format(total_score, right_score))
        logging.info("right match {}({})".format(right_match, right_match/float(len(indexes))))


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    gpath1 = r"./multiplex_embeddings_data/IPv4_IPv6 Internet/l_1-links.txt"
    gpath2 = r"./multiplex_embeddings_data/IPv4_IPv6 Internet/l_2-links.txt"
    gpath1 = r"./multiplex_embeddings_data/Drosophila Melanogaster/l_1-links.txt"
    gpath2 = r"./multiplex_embeddings_data/Drosophila Melanogaster/l_2-links.txt"
    gpath1 = r"./multiplex_embeddings_data/Air_Train/l_1.txt"
    gpath2 = r"./multiplex_embeddings_data/Air_Train/l_2.txt"
    v = ShortestPathLengthVector(gpath1, gpath2)
    v.km_match()