import networkx as nx
from . import utils

class MigrationGraph:

    def __init__(self, raw):
        G_raw = raw['aug_mig_graph']
        self.graph = nx.MultiDiGraph()
        for i in G_raw.values():
            for s, t in i:
                if s != t:
                    self.graph.add_edge(s, t)
        self.sites = [i for i in self.graph.nodes]
        self.n_mig = raw['n_mig']
        self.n_comig = raw['n_comig']
        if 'n_seedingsites' in raw:
            self.n_seedingsites = raw['n_seedingsites']

    def has_migration(self, a, b):
        return self.graph.number_of_edges(a,b) > 0
    
    def n_migrations(self, a, b):
        return self.graph.number_of_edges(a,b)
    
    def migration_edges(self):
        G_digraph = nx.DiGraph(self.graph)
        for u, v in G_digraph.edges:
            yield u, v
    
    def migration_pattern(self):
        G_digraph = nx.DiGraph(self.graph)
        clonality = 'm'
        for edge in G_digraph.edges:
            if self.graph.n_migrations(edge[0], edge[1]) > 1:
                clonality = 'p'
                break
        if nx.is_tree(G_digraph):
            return clonality, 'S'
        elif nx.is_directed_acyclic_graph(G_digraph):
            return clonality, 'M'
        else:
            return clonality, 'R'
    
    def write_graph(self, filename):
        G_digraph = nx.DiGraph(self.graph)
        with open(filename, 'w+') as f:
            for edge in G_digraph.edges:
                f.write(f'{edge[0]} {edge[1]} {self.graph.number_of_edges(edge[0], edge[1])}\n')

    def draw(self, colormap=None, colormap_file=None):
        import graphviz as gv
        if colormap is None:
            if colormap_file is None:
                colormap = utils.get_colormap(self.sites)
            else:
                colormap = utils.process_colormap_file(colormap_file)
        g = gv.Digraph(node_attr={'shape': 'box', 'penwidth': '3', 'colorscheme': 'set19'}, edge_attr={'penwidth': '3', 'colorscheme': 'set19'})
        for s in self.graph.nodes:
            g.node(s, color=str(colormap[s]))
        for s, t, _ in self.graph.edges:
            g.edge(s, t, color=f"{colormap[s]};0.5:{colormap[t]}")
        return g

    def write_dot(self, filename, colormap=None, colormap_file=None):
        if colormap is None:
            if colormap_file is None:
                colormap = utils.get_colormap(self.sites)
            else:
                colormap = utils.process_colormap_file(colormap_file)
        node_index = {j:i for i,j in enumerate(self.graph.nodes)}
        with open(filename, 'w+') as f:
            f.write('digraph G {\n')
            for s in self.graph.nodes:
                f.write(f'\t{node_index[s]} [shape=box,penwidth=3,colorscheme=set19,color={colormap[s]},' + f'label="{s}"]\n')
            for s, t, _ in self.graph.edges:
                f.write(f'\t{node_index[s]} -> {node_index[t]} [penwidth=3,colorscheme=set19,' +
                            f'color="{colormap[s]};0.5:{colormap[t]}"]\n')
            f.write('}\n')
