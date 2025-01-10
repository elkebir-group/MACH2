import networkx as nx
from collections import Counter
from . import utils

class MigrationGraph:

    def __init__(self, refinement):
        self.refinement = refinement
        self._graph = nx.MultiDiGraph()
        for u, v in refinement.edges:
            if refinement.get_label(u) != refinement.get_label(v):
                self._graph.add_edge(refinement.get_label(u), refinement.get_label(v))

    def has_migration(self, a, b):
        return self._graph.number_of_edges(a,b) > 0
    
    def n_migrations(self, a, b):
        return self._graph.number_of_edges(a,b)
    
    def get_each_migration(self):
        return Counter([(s, t) for s, t, _ in self._graph.edges])
    
    def get_migrations(self):
        G_digraph = nx.DiGraph(self._graph)
        for u, v in G_digraph.edges:
            yield u, v

    def get_seeding_locations(self):
        return set([u for u, _ in self.get_migrations()])
    
    def migration_pattern(self):
        G_digraph = nx.DiGraph(self._graph)
        clonality = 'm'
        for edge in G_digraph.edges:
            if self.n_migrations(edge[0], edge[1]) > 1:
                clonality = 'p'
                break
        if nx.is_tree(G_digraph):
            return clonality, 'S'
        elif nx.is_directed_acyclic_graph(G_digraph):
            return 'p', 'M'
        else:
            return 'p', 'R'
    
    def write_graph(self, filename):
        G_digraph = nx.DiGraph(self._graph)
        with open(filename, 'w') as f:
            for edge in G_digraph.edges:
                f.write(f'{edge[0]}\t{edge[1]}\t{self._graph.number_of_edges(edge[0], edge[1])}\n')

    def draw(self):
        import graphviz as gv
        colormap = self.refinement.unrefined_tree._colormap
        g = gv.Digraph(node_attr={'shape': 'box', 'penwidth': '3', 'colorscheme': 'set19'}, edge_attr={'penwidth': '3', 'colorscheme': 'set19'})
        for s in self._graph.nodes:
            g.node(s, color=str(colormap[s]))
        for s, t, _ in self._graph.edges:
            if s != t:
                g.edge(s, t, color=f"{colormap[s]};0.5:{colormap[t]}")
        return g
    
    def _repr_mimebundle_(self, include=None, exclude=None, **_):
        return self.draw()._repr_mimebundle_(include=include, exclude=exclude)