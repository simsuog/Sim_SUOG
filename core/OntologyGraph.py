#
#  Copyright (c) 2023 - The SUOG Project.
#  @Author: Mirna El Ghosh
#

from abc import abstractmethod, ABCMeta
import networkx as nx
from core.SUOGOntology import SUOGOnto
import numpy as np



class OntologyGraph:

    __metaclass__ = ABCMeta

    @abstractmethod
    def define_graph(self):
        """return nodes, labels, and hierarchical edges"""
        pass


class SUOGOntologyGraph(OntologyGraph):

    """To build the hierarchical DAG of the SUOG ontology (owl file)"""

    def __init__(self, src):
        self._ontology = SUOGOnto(src)
        self._root = self._ontology._root

    def define_graph(self):
        nodes = {n for n in self._ontology._listclass}
        labels = tuple([self._ontology.labelClassOf(value) for i, value in enumerate(self._ontology._listclass)])
        edges = []
        for i, node in enumerate(nodes):
            children = self._ontology.subClassesOf(node)
            children = [child for child in children if child in nodes]
            for child in children:
               edges.append((node, child))

        return nodes, labels, edges


class Taxonomy(object):

    def __init__(self, onto):
        self._nodes, self._labels, self._edges = onto.define_graph()
        self._node2id = {value: i for i, value in enumerate(self._nodes)}
        self._root = onto._root
        self._taxonomy = nx.DiGraph()
        self._hyponyms = {}
        self._hypernyms = {}
        self._leavesMax = {}
        self._children = {}
        self.build_graph()

    def build_graph(self):
        parents, children = zip(*self._edges)
        parents_set = set(parents)
        children_set = set(children)
        self._children = children_set
        not_children = []
        not_parent = []

        #parents that are not in children set

        for p in parents_set:
            if p not in children_set:
                not_children.append(p)

        # children that are not in parent set
        for p in children_set:
            if p not in parents_set:
                not_parent.append(p)
        self._leavesMax = not_parent

        self._taxonomy.add_nodes_from(self._nodes)

        # add taxonomical edges
        for parent, child in self._edges:
            self._taxonomy.add_edge(parent, child)

        # hyponyms and hypernyms
        for parent, child in self._edges:
            self._hyponyms.setdefault(parent,[]).append(child)

        for parent, child in self._edges:
            self._hypernyms.setdefault(child, []).append(parent)

    def root_children(self):
        return self.hyponyms(self._root)

    def is_directed(self):
        return self._taxonomy.is_directed()

    def max_leaves(self):
        return self._leavesMax.__len__()

    def hyponyms(self, node):
        return self._hyponyms[node] if node in self._hyponyms else []

    def hypernyms(self, node):
        return self._hypernyms[node] if node in self._hypernyms else []

    def is_parent(self, c1, c2):

        return c1 in self.hypernyms(c2)

    def shortest_path_length(self, node1, node2):
        length = 0
        if nx.has_path(self._taxonomy, node1, node2):
            length = (nx.shortest_path_length(self._taxonomy, node1, node2))
        return length

    def depth(self, node):
        print('depth of ', node,': ', len(self.shortest_path_length(self._root, node)))

        return len(self.shortest_path_length(self._root, node))

    def max_depth(self, onto):
        depths = []
        for c in onto.nodes:
           # print(c,self._root,self.shortest_path_length(c, self._root))
            depths.append(self.shortest_path_length(c, self._root))
        return max(depths)

    def lca(self, c1, c2):
        l = None
        if nx.lowest_common_ancestor is not None:
            l = nx.lowest_common_ancestor(self._taxonomy, c1, c2)
        return l

    def leaves(self, c):

        return [h for h in self.hyponyms(c) if h in self._leavesMax]

    def leaves_all(self, c):

        return [h for h in self.descendants(c) if h in self._leavesMax]

    def IC(self, x):

        """IC is computed using Sanchez's equation"""

        ic = -np.log((((self.leaves_all(x).__len__()) / (self.ancestors(x).__len__() + 1)) + 1) / (self.max_leaves() + 1))

        return (ic)

    def descendants(self, c):
        return nx.descendants(self._taxonomy, c)

    def ancestors(self, c):
        return nx.ancestors(self._taxonomy, c)

    def common_ancestors(self, c1, c2):
        ac1 = nx.ancestors(self._taxonomy,c1)
        ac2 = nx.ancestors(self._taxonomy,c2)
        listca = [a for a in ac1 if a in ac2]

        return listca

    def siblings(self, c):
        sib = []
        for parent in self.hypernyms(c):
            for child in self.hyponyms(parent):
                if child != c:
                    sib.append(child)
        return sib

    def mica(self, c1, c2):
        listca = self.common_ancestors(c1, c2)
        listic = []
        listic2 = []
        for l in listca:
            ic = self.IC(l)
            listic.append(ic)
            listic2.append(l)
        maxic = max(listic)
        index = listic.index(maxic)
        mica = listic2[index]

        return mica

    def sim_resnik(self, c1, c2):

        sim = self.IC(self.mica(c1, c2))

        return np.round(sim, 2)

    def sim_lin(self, c1, c2):

        sim = np.round(((2 * self.sim_resnik(c1, c2)) / (self.IC(c1) + self.IC(c2))), 2)

        return sim

    def sim_wu_palmer(self, c1 , c2):

        sim = np.round(((2 * self.depth(self.lca(c1, c2))) / (self.depth(c1) + self.depth(c2))), 2)

        return sim

    def sim_ic(self, c1, c2):

        sim = self.sim_lin(c1, c2) * (1 - (1 / (1 + self.sim_resnik(c1, c2))))
        return sim

    def ca(self, c1, c2):
        return self.common_ancestors(c1,c2)

    def sim_suog(self, c1, c2):
        if c1 == c2:
          sim = 1
        else:
            if self.lca(c1, c2) is None:
                sim = 0
            else:
                dist1 = self.shortest_path_length(self.lca(c1, c2), c1)
                dist2 = self.shortest_path_length(self.lca(c1, c2), c2)
                dist = (dist1+dist2)/2
                ic = self.IC(self.lca(c1, c2))
                if dist > ic:
                    sim = 0
                else:
                    sim = 1-np.round(dist/ic, 3)
        return sim

    def sim_suog_t(self, c1, c2):

        if c1 == c2:
            sim = 1
        else:
            dist1 = self.shortest_path_length(self.lca(c1, c2), c1)
            dist2 = self.shortest_path_length(self.lca(c1, c2), c2)
            dist = (dist1+dist2)/2
            ic = self.IC(self.lca(c1, c2))
            sim = np.abs(1- np.round((2*(1+dist)) / ic, 3))

        return sim