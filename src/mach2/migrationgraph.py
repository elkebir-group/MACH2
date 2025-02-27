import networkx as nx
from collections import Counter
import graphviz as gv

class MigrationGraph:
    """
    A class representing a migration graph derived from a Refinement object.
    This graph captures the migration patterns between locations in the refinement.
    """
    def __init__(self, refinement):
        """
        Initializes the MigrationGraph with a Refinement object.
        
        :param refinement: A Refinement object representing the refined tree.
        """
        self.refinement = refinement
        self._graph = nx.MultiDiGraph()
        for u, v in refinement.edges:
            if refinement.get_label(u) != refinement.get_label(v):
                self._graph.add_edge(refinement.get_label(u), refinement.get_label(v))

    def has_migration(self, a, b):
        """
        Checks if there is a migration from location `a` to location `b`.
        
        :param a: The source location.
        :param b: The target location.
        :return: True if there is at least one migration, False otherwise.
        """
        return self._graph.number_of_edges(a, b) > 0
    
    def n_migrations(self, a, b):
        """
        Returns the number of migrations from location `a` to location `b`.
        
        :param a: The source location.
        :param b: The target location.
        :return: The number of migrations as an integer.
        """
        return self._graph.number_of_edges(a, b)
    
    def get_each_migration(self):
        """
        Returns a Counter object with the count of each unique migration edge.
        
        :return: A Counter object where keys are tuples (source, target) and values are counts.
        """
        return Counter([(s, t) for s, t, _ in self._graph.edges])
    
    def get_migrations(self):
        """
        Yields all unique migration edges in the graph.
        
        :return: A generator yielding tuples (source, target) representing migrations.
        """
        G_digraph = nx.DiGraph(self._graph)
        for u, v in G_digraph.edges:
            yield u, v

    def get_seeding_locations(self):
        """
        Returns the set of locations that act as sources for migrations.
        
        :return: A set of seeding locations.
        """
        return set([u for u, _ in self.get_migrations()])
    
    def migration_pattern(self):
        """
        Determines the migration pattern of the graph.
        
        :return: A tuple (clonality, pattern) where:
                 - clonality is 'm' (monoclonal) or 'p' (polyclonal).
                 - pattern is 'S' (single-source), 'M' (multi-source), or 'R' (reticulate).
        """
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
        """
        Writes the migration graph to a file.
        
        :param filename: The path to the output file.
        """
        G_digraph = nx.DiGraph(self._graph)
        with open(filename, 'w') as f:
            for edge in G_digraph.edges:
                f.write(f'{edge[0]}\t{edge[1]}\t{self._graph.number_of_edges(edge[0], edge[1])}\n')

    def draw(self, filename=None):
        """
        Draws the migration graph using Graphviz.
        
        :param filename: Optional path to save the rendered graph.
        :return: A Graphviz Digraph object.
        """
        colormap = self.refinement.unrefined_tree._colormap
        g = gv.Digraph(node_attr={'shape': 'box', 'penwidth': '3', 'colorscheme': 'set19'}, edge_attr={'penwidth': '3', 'colorscheme': 'set19'})
        for s in self._graph.nodes:
            g.node(s, color=str(colormap[s]))
        for s, t, _ in self._graph.edges:
            if s != t:
                g.edge(s, t, color=f"{colormap[s]};0.5:{colormap[t]}")
        if filename is not None:
            g.render(filename)
        return g
    
    def _repr_mimebundle_(self, include=None, exclude=None, **_):
        """
        Returns the MIME bundle representation of the graph for visualization.
        
        :param include: Optional set of MIME types to include.
        :param exclude: Optional set of MIME types to exclude.
        :return: A dictionary containing the MIME bundle.
        """
        return self.draw()._repr_mimebundle_(include=include, exclude=exclude)