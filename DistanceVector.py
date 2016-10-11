# -*- coding: utf-8 -*-
import networkx as nx
import numpy as np
import random
import itertools
import logging
from collections import namedtuple
from math import pi, fabs, exp, acosh, cosh, sinh, cos, log


HybolicPara = namedtuple("HybolicPara", ["radius", "theta"])


def construct_hypermap_from_file(path):
    """read HyperMap program output."""
    r = {}
    with open(path, "r") as f:
        for line in f.readlines():
            id, theta, radius = line.split()
            p = HybolicPara(radius=float(radius), theta=float(theta))
            r[id] = p
    logging.debug("{} nodes and hyperbolic parameter loaded from {}.".format
                (len(r), path))
    return r


def vector_distance(x, y):
    x, y = np.array(x), np.array(y)
    lx, ly = np.sqrt(x.dot(x)), np.sqrt(y.dot(y))
    cos_angle = x.dot(y) / (lx * ly)
    return np.arccos(cos_angle)


class HybolicDistanceVector(object):

    def __init__(self, layer_path, hybolic_path, hybolic_T_R):
        assert(len(layer_path) == len(hybolic_path) == len(hybolic_T_R) >= 2)
        self.layer = range(len(layer_path))
        self.graphs = []
        self.hyperspace = []
        self.hyperspace_TR = hybolic_T_R
        self.revealed_nodes = {}
        self.same_nodes = {}
        for l in self.layer:
            self.graphs.append(nx.read_edgelist(layer_path[l]))
            logging.debug("{} read {}.".format(l, layer_path[l]))
            self.hyperspace.append(construct_hypermap_from_file(hybolic_path[l]))
            logging.debug("{} read {}.".format(l, hybolic_path[l]))

        for i, j in itertools.combinations(self.layer, 2):
            self.random_select_revealed_node(0.9, i, j)

    def hybolic_distance(self, index, x, y):
        """hybolic space distance between point x and y."""
        r1, t1 = self.hyperspace[index][x]
        r2, t2 = self.hyperspace[index][y]
        delta_t = (pi - fabs(pi - fabs(t1 - t2))) + 0.0001
        d = r1 + r2 + 2.0 * log(delta_t / 2.0)
        #d = acosh(cosh(r1) * cosh(r2) - sinh(r1) * sinh(r2) * cos(delta_t))
        T, R = self.hyperspace_TR[index]
        #p = 100.0 / (1.0 + exp((d - R) / (2.0 * T)))  # percent
        p = 1.0 if d >= R else 0.0
        return p

    def random_select_revealed_node(self, alpha, index1, index2):
        """select some nodes as revealed nodes randomly,
        alpha is percent of all same nodes."""
        same_nodes = set(self.graphs[index1].nodes()) & set(self.graphs[index2].nodes())
        s = int(alpha * len(same_nodes))
        logging.info("graph {}-{} random revealed nodes {}/{}.".format
                     (index1, index2, s, len(same_nodes)))
        revealed_nodes = random.sample(same_nodes, s)
        self.revealed_nodes[(index1, index2)] = set(revealed_nodes)
        self.same_nodes[(index1, index2)] = same_nodes
        return revealed_nodes

    def unmatched_node_distance_vector(self, index1, index2):
        unmatched_nodes = self.same_nodes[(index1, index2)] - self.revealed_nodes[(index1, index2)]
        dv1, dv2 = {}, {}
        for n in unmatched_nodes:
            r1, r2 = [], []
            for r in self.revealed_nodes[(index1, index2)]:
                r1.append(self.hybolic_distance(index1, n, r))
                r2.append(self.hybolic_distance(index2, n, r))
            dv1[n] = r1
            dv2[n] = r2
        return dv1, dv2

    def cal_match_node(self, index1, index2):
        dv1, dv2 = self.unmatched_node_distance_vector(index1, index2)
        for n1, v1 in dv1.iteritems():
            r = {}
            for n2, v2 in dv2.iteritems():
                d12 = vector_distance(v1, v2)
                r[n2] = d12
                if n2 == n1:
                    logging.debug("{} right match score {}.".format(n1, d12))
            sr = sorted(r.iteritems(), key=lambda k: k[1], reverse=False)
            logging.info("{} 1-mostly possible match {} with score {}.".format
                        (n1, sr[0][0], sr[0][1]))
            for i, (n_, s_) in enumerate(sr):
                if n_ == n1:
                    logging.info("{} right match is at {} pos with score {}.".format
                        (n1, i, s_))


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    gpath1 = r"./multiplex_embeddings_data/Drosophila Melanogaster/l_1-links.txt"
    gpath2 = r"./multiplex_embeddings_data/Drosophila Melanogaster/l_2-links.txt"
    hpath1 = r"./multiplex_embeddings_data/Drosophila Melanogaster/l_1-coordinates.txt"
    hpath2 = r"./multiplex_embeddings_data/Drosophila Melanogaster/l_2-coordinates.txt"
    c = HybolicDistanceVector([gpath1, gpath2], [hpath1, hpath2], [(0.1, 12.4915), (0.1, 12.5384)])
    c.cal_match_node(0, 1)

