# -*- coding: utf-8 -*-
from collections import namedtuple, defaultdict
import itertools
import logging
import random
from math import pi, log, fabs, exp
import matplotlib.pyplot as plt
import networkx as nx


HybolicPara = namedtuple("HybolicPara", ["radius", "theta"])



class HyperMapAnalysis(object):

    def __init__(self, merge_g_path, hypermap_path, T):
        """multi_g is MergeMultiGraph"""
        self.g = nx.read_edgelist(merge_g_path)
        self.R = 2.0 * log(self.g.number_of_nodes())
        logging.info("{} nodes and {} edges loaded from {}, R={}.".format
                    (self.g.number_of_nodes(), self.g.number_of_edges(), merge_g_path, self.R))
        self.T = T
        self.hypermap = self.construct_hypermap_from_file(hypermap_path)

    def construct_hypermap_from_file(self, path):
        """read HyperMap program output."""
        r = {}
        with open(path, "r") as f:
            for line in f.readlines():
                id, theta, radius = line.split()
                p = HybolicPara(radius=float(radius), theta=float(theta))
                r[id] = p
        logging.info("{} nodes and hyperbolic parameter loaded from {}.".format
                    (len(r), path))
        return r

    def connected_nodes_correlation(self):
        """show connected nodes correlation."""
        connected_nodes = self.g.edges()
        r = []
        for x, y in connected_nodes:
            d = self.hybolic_distance(x, y)
            r.append(d)
        return r

    def unconnected_nodes_correlation(self):
        """show unconnected nodes correlation."""
        n = self.g.number_of_edges()
        random_edges = random.sample(list(itertools.combinations(self.g.nodes(), 2)), n)
        r = []
        for i, (x, y) in enumerate(random_edges):
            if self.g.has_edge(x, y):
                continue
            d = self.hybolic_distance(x, y)
            r.append(d)
            if i >= n:
                break
        logging.info("{} random unconnected nodes.".format(i))
        return r

    def unmatched_nodes_correlation(self):
        """matched nodes except revealed ones."""
        unmatched_nodes = []
        for x in self.g.nodes():
            _, rx = x.split("_")
            for y in self.g.nodes():
                if y <= x:
                    continue
                _, ry = y.split("_")
                if rx == ry:
                    unmatched_nodes.append((x, y))
        logging.info("{} unmatched nodes.".format
                    (len(unmatched_nodes)))
        r = []
        for x, y in unmatched_nodes:
            d = self.hybolic_distance(x, y)
            r.append(d)
        return r

    def hybolic_distance(self, x, y):
        """hybolic space distance between point x and y."""
        r1, t1 = self.hypermap[x]
        r2, t2 = self.hypermap[y]
        delta_t = (pi - fabs(pi - fabs(t1 - t2))) + 0.0001
        logging.debug("{}-{} hybolic delta theta: {}.".format
                     (x, y, delta_t))
        d = r1 + r2 + 2.0 * log(delta_t / 2.0)
        p = 100.0 / (1.0 + exp((d - self.R) / (2.0 * self.T)))  # percent
        return d, p

    def list_statisic(self, l, step=1):
        """l is a list with number, return it's statisic by step."""
        rd, rp = defaultdict(int), defaultdict(int)
        for d, p in l:
            x = int(d / step) * step
            rd[x] += 1
            x = int(p / step) * step
            rp[x] += 1
        return rd, rp


def plot_bar(d, subp, title="", xlabel="", ylabel=""):
    p = plt.subplot(subp)
    p.bar(d.keys(), d.values())
    p.set_title(title)
    p.set_xlabel(xlabel)
    p.set_ylabel(ylabel)
    return p


def plot2():
    path = "arXiv/l_8"
    a = HyperMapAnalysis("./multiplex_embeddings_data/{}.txt".format(path),
    "./multiplex_embeddings_data/{}_coordinates.txt".format(path), T=0.05)
    connected_dp = a.connected_nodes_correlation()
    unconnected_dp = a.unconnected_nodes_correlation()

    connected_d, connected_p = a.list_statisic(connected_dp)
    unconnected_d, unconnected_p = a.list_statisic(unconnected_dp)
    plot_bar(connected_d, 221, "connected d", "dis", "num")
    plot_bar(connected_p, 222, "connected p", "p", "num")
    plot_bar(unconnected_d, 223, "unconnected d", "dis", "num")
    plot_bar(unconnected_p, 224, "unconnected p", "p", "num")
    plt.savefig("{}_{}.png".format(path.replace("/", "_"), int(a.R)))
    plt.show()


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    plot2()





