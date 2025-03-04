import networkx as nx
import graphviz as gv
from .migrationgraph import MigrationGraph
from . import utils
from collections import defaultdict


class Tree:
    """
    A class representing a tree structure.

    Class Variables:
        _tree (nx.DiGraph): A directed graph representing the tree structure.
        leaves (list): A list of leaf nodes in the tree.
        root (str): The root node of the tree.
        nodes (list): A list of all nodes in the tree.
        edges (list): A list of all edges in the tree.
        paths (dict): A dictionary mapping each node to its path from the root.
        max_height (int): The maximum height (longest path) of the tree.
    """
    def __init__(self, edges):
        """
        Initializes the tree from a given list of edges.

        :param edges: List of tuples representing the edges of the tree.
        :raises ValueError: If the provided edges do not form a valid tree.
        """
        self._tree = nx.DiGraph()
        for edge in edges:
            self._tree.add_edge(edge[0], edge[1])
        if not nx.is_tree(self._tree):
            raise ValueError('Not a valid tree')

        self.leaves = [i for i in self._tree.nodes if self._tree.out_degree(i) == 0]
        self.root = [i for i in self._tree.nodes if self._tree.in_degree(i) == 0][0]
        self.nodes = [v for v in self._tree.nodes]
        self.edges = [e for e in self._tree.edges]

        self.paths = {}
        self.max_height = 0
        for u in self.nodes:
            self.paths[u] = []
            path_u = nx.shortest_path(self._tree, self.root, u)
            for ii in range(len(path_u) - 1):
                self.paths[u].append((path_u[ii], path_u[ii+1]))
            if len(self.paths[u]) > self.max_height:
                self.max_height = len(self.paths[u])

    def __eq__(self, other):
        """
        Checks if two trees are equal based on their edges.

        :param other: Another Tree object to compare with.
        :return: True if the trees are equal, False otherwise.
        """
        return set(self.edges) == set(other.edges)
    
    def __ne__(self, other):
        """
        Checks if two trees are not equal based on their edges.

        :param other: Another Tree object to compare with.
        :return: True if the trees are not equal, False otherwise.
        """
        return set(self.edges) != set(other.edges)
    
    def __hash__(self):
        """
        Returns a hash value for the tree.

        :return: A hash value based on the tree's edges.
        """
        return hash(frozenset(self.edges))

    def from_file(filename):
        """
        Creates a Tree instance from a file containing edges.

        :param filename: Path to the file containing edge data.
        :return: A Tree instance.
        :raises ValueError: If the file format is incorrect.
        """
        with open(filename, 'r') as f:
            edges = []
            for edge in f.readlines():
                try:
                    a, b = edge.split()
                except:
                    raise ValueError('Ill-formatted input tree file')
                edges.append((a, b))
        return Tree(edges)
    
    def get_parent(self, node):
        """
        Returns the parent of a given node (string).

        :param node: The node (string) whose parent is to be retrieved.
        :return: The parent node (string), or None if the node is the root.
        """
        try:
            return [k for k in self._tree.in_edges(nbunch=node)][0][0]
        except:
            return None
    
    def get_children(self, node):
        """
        Returns the children of a given node (string).

        :param node: The node (string) whose children are to be retrieved.
        :return: A list of child nodes (string).
        """
        return [k[1] for k in self._tree.out_edges(nbunch=node)]
    
    def get_descendant_leaves(self, node):
        """
        Returns all descendant leaves of a given node (string).

        :param node: The node (string) whose descendant leaves are to be retrieved.
        :return: A list of descendant leaf nodes (string).
        """
        if node in self.leaves:
            return [node]
        desc = []
        for v in self.get_children(node):
            desc += self.get_descendant_leaves(v)
        return desc
    
    def get_descendants(self, node):
        """
        Returns all descendants of a given node (string).

        :param node: The node (string) whose descendants are to be retrieved.
        :return: A list of descendant nodes (string).
        """
        desc = []
        for v in self.get_children(node):
            desc += [v] + self.get_descendants(v)
        return desc
    
    def get_path(self, node1, node2):
        """
        Returns the shortest path between two nodes.

        :param node1: The starting node (string).
        :param node2: The target node (string).
        :return: A list of edges (tuple of strings) representing the shortest path, or None if no path exists.
        """
        try:
            return nx.shortest_path(self._tree, node1, node2)
        except:
            return None

    def write_tree(self, filename):
        """
        Writes the tree edges to a file.

        :param filename: Path to the output file.
        """
        with open(filename, 'w') as f:
            for edge in self.edges:
                f.write(f'{edge[0]}\t{edge[1]}\n')

    def draw(self):
        """
        Draws the tree using Graphviz.

        :return: A Graphviz Digraph object representing the tree.
        """
        t = gv.Digraph(node_attr={'penwidth': '3'}, edge_attr={'penwidth': '3'})
        t_leaves = gv.Digraph(graph_attr={'rank': 'same'}, node_attr={'penwidth': '3'})
        for l in self.leaves:
            t_leaves.node(l, f"{l}")
        t.subgraph(t_leaves)
        for i, j in self.edges:
            t.edge(i, j)
        return t
    
    def _repr_mimebundle_(self, include=None, exclude=None, **_):
        """
        Private method to return the MIME bundle representation of the tree for visualization.
        """
        return self.draw()._repr_mimebundle_(include=include, exclude=exclude)
    

class MultiLabeledTree(Tree):
    """
    A class representing a tree with multiple labels assigned to nodes.

    Class Variables:
        _colormap (dict): A dictionary mapping labels to colors (represented by index in Brewer color scheme Set19).
        locations (set): A set of all unique locations (labels) in the tree.
    """
    def __init__(self, edges, labels, colormap=None):
        """
        Initializes the MultiLabeledTree with edges, labels, and an optional colormap.

        :param edges: List of tuples of strings representing the edges of the tree.
        :param labels: Dictionary mapping nodes (strings) to their labels (string).
        :param colormap: Optional dictionary mapping labels (string) to colors (integer).
        """
        super().__init__(edges)
        attrs = {}
        for node in self.nodes:
            if node in labels:
                attrs[node] = {'label': set(labels[node])}
            else:
                attrs[node] = {'label': set()}
        nx.set_node_attributes(self._tree, attrs)
        self.locations = set(ss for s in labels for ss in labels[s])

        if colormap is None:
            self._colormap = utils.get_colormap(self.locations)
        else:
            self._colormap = colormap

    def __eq__(self, other):
        """
        Checks if two MultiLabeledTrees are equal based on edges and labels.

        :param other: Another MultiLabeledTree object to compare with.
        :return: True if the trees are equal, False otherwise.
        """
        return super().__eq__(other) and set([(u, frozenset(self.get_labels(u))) for u in self.nodes]) == set([(u, frozenset(other.get_labels(u))) for u in other.nodes])
    
    def __ne__(self, other):
        """
        Checks if two MultiLabeledTrees are not equal based on edges and labels.

        :param other: Another MultiLabeledTree object to compare with.
        :return: True if the trees are not equal, False otherwise.
        """
        return not self == other
    
    def __hash__(self):
        """
        Returns a hash value for the MultiLabeledTree.

        :return: A hash value based on the tree's edges and labels.
        """
        return hash((frozenset(self.edges), frozenset([(u, frozenset(self.get_labels(u))) for u in self.nodes])))

    def from_files(treefile, labelfile, colormap_file=None):
        """
        Creates a MultiLabeledTree instance from files containing edges and labels.

        :param treefile: Path to the file containing edge data.
        :param labelfile: Path to the file containing label data.
        :param colormap_file: Optional path to the file containing colormap data.
        :return: A MultiLabeledTree instance.
        :raises ValueError: If the file format is incorrect.
        """
        with open(treefile, 'r') as f:
            edges = []
            for edge in f:
                try:
                    a, b = edge.split()
                except:
                    raise ValueError('Ill-formatted input tree file')
                edges.append((a, b))
        with open(labelfile, 'r') as f:
            labels = {}
            for line in f:
                l = line.split()
                labels[l[0]] = l[1:]
        colormap = None
        if colormap_file is not None:
            colormap = utils.process_colormap_file(colormap_file)
        return MultiLabeledTree(edges, labels, colormap)
    
    def get_labels(self, node):
        """
        Returns the labels of a given node.

        :param node: The node (string) whose labels are to be retrieved.
        :return: A set of labels (string) associated with the node.
        """
        return self._tree.nodes[node]['label']
    
    def add_label(self, node, label):
        """
        Adds a label to a given node.

        :param node: The node (string) to which the label is to be added.
        :param label: The label (string) to add.
        """
        if label not in self.locations:
            self.locations.add(label)
        self._tree.nodes[node]['label'].add(label)

    def write_labelings(self, filename):
        """
        Writes the node labels to a file.

        :param filename: Path to the output file.
        """
        with open(filename, 'w') as f:
            for node in self.nodes:
                labels = self.get_labels(node)
                if len(labels) > 0:
                    f.write(f'{node}')
                    for label in labels:
                        f.write(f'\t{label}')
                    f.write('\n')

    def set_colormap(self, colormap):
        """
        Sets the colormap for the tree.

        :param colormap: A dictionary mapping labels (string) to colors (integer).
        """
        self._colormap = colormap

    def set_time_of_observation(self, time_of_observation):
        """
        Sets the time of observation for the tree.

        :param time_of_observation: The time of observation to set.
        """
        self.time_of_observation = time_of_observation

    def draw(self):
        """
        Draws the MultiLabeledTree using Graphviz.

        :return: A Graphviz Digraph object representing the tree.
        """
        t = gv.Digraph(node_attr={'penwidth': '3', 'colorscheme': 'set19'}, edge_attr={
                       'penwidth': '3', 'colorscheme': 'set19'})
        t_leaves = gv.Digraph(graph_attr={'rank': 'same'}, node_attr={
                              'penwidth': '3', 'colorscheme': 'set19'})
        for l in self.leaves:
            t_leaves.node(l, f"{l}\n{self.get_labels(l)}", color=str(0))
        t.subgraph(t_leaves)
        for v in self.nodes:
            if v not in self.leaves:
                t.node(v, f"{v}\n{self.get_labels(v)}", color=str(0))
        for i, j in self.edges:
            t.edge( i, j, color=str(0))
        return t


class Refinement(Tree):
    """
    A class representing a refined tree with additional attributes.

    Class Variables:
        unrefined_tree (MultiLabeledTree): The original unrefined tree.
        node_of_origin (dict): A dictionary mapping nodes to their origin nodes.
        unobserved_clones (list): A list of unobserved clones in the refinement.
        migrations (list): A list of migration edges in the refinement.
        comigrations (defaultdict): A dictionary mapping timestamps to comigration events.
        seeding_locations (list): A list of seeding locations in the refinement.
        timestamps (dict): A dictionary mapping edges to their timestamps.
    """
    def __init__(self, edges, labels, unrefined_tree, node_of_origin=None, timestamps=None):
        """
        Initializes the Refinement tree with edges, labels, and additional attributes.

        :param edges: List of tuples representing the edges of the tree.
        :param labels: Dictionary mapping nodes to their labels.
        :param unrefined_tree: The original unrefined tree.
        :param node_of_origin: Dictionary mapping nodes to their origin nodes.
        :param timestamps: Dictionary mapping edges to their timestamps.
        """
        super().__init__(edges)
        attrs = {i: {'label': j} for (i, j) in labels.items()}
        nx.set_node_attributes(self._tree, attrs)
        self.locations = set(labels.values())
        self.unrefined_tree = unrefined_tree
        if node_of_origin is None:
            node_of_origin = {u : u.split('^')[0] for u in self.nodes}
        self.node_of_origin = node_of_origin
        self.unobserved_clones = [up for up in self.nodes if self.get_label(up) not in unrefined_tree.get_labels(node_of_origin[up])]
        self.migrations = [(u, v) for (u, v) in self.edges if self.get_label(u) != self.get_label(v)]
        if timestamps is None:
            try:
                self._get_comigrations()
            except ModuleNotFoundError:
                print('Gurobi required for inferring temporally consistent comigrations not found. \
                      Comigrations inferred by greedy approach may be temporally inconsistent.')
                self._greedy_comigrations()
        else:
            self.timestamps = timestamps
            self.comigrations = defaultdict(list)
            for uv, t in timestamps.items():
                self.comigrations[t].append(uv)
        self.seeding_locations = list(set([attrs[u]['label'] for (u, _) in self.migrations]))

    def __eq__(self, other):
        """
        Checks if two Refinement trees are equal based on edges and labels.

        :param other: Another Refinement object to compare with.
        :return: True if the trees are equal, False otherwise.
        """
        return super().__eq__(other) and set([(u, self.get_label(u)) for u in self.nodes]) == set([(u, other.get_label(u)) for u in other.nodes])
    
    def __ne__(self, other):
        """
        Checks if two Refinement trees are not equal based on edges and labels.

        :param other: Another Refinement object to compare with.
        :return: True if the trees are not equal, False otherwise.
        """
        return not self == other
    
    def __hash__(self):
        """
        Returns a hash value for the Refinement tree.

        :return: A hash value based on the tree's edges and labels.
        """
        return hash((frozenset(self.edges), frozenset([(u, self.get_label(u)) for u in self.nodes])))

    def from_refinement_graph(unrefined_tree, refinement_graph):
        """
        Creates a Refinement instance from an unrefined tree and a refinement graph.

        :param unrefined_tree: The original unrefined tree.
        :param refinement_graph: The refinement graph.
        :return: A Refinement instance.
        """
        edges = []
        labels = {}
        node_of_origin = {}
        for v in unrefined_tree.nodes:
            if v in refinement_graph:
                for s, t in refinement_graph[v]:
                    labels[f'{v}^{s}'] = s
                    labels[f'{v}^{t}'] = t
                    node_of_origin[f'{v}^{s}'] = v
                    node_of_origin[f'{v}^{t}'] = v
                    edges.append((f'{v}^{s}', f'{v}^{t}'))
            else:
                node_of_origin[v] = v
        for u, v in unrefined_tree.edges:
            s, t = refinement_graph[(u,v)][0]
            up = f'{u}^{s}' if f'{u}^{s}' in node_of_origin else u
            vp = f'{v}^{t}' if f'{v}^{t}' in node_of_origin else v
            edges.append((up, vp))
            labels[up] = s
            labels[vp] = t
        return Refinement(edges, labels, unrefined_tree, node_of_origin)
    
    def from_files(treefile, labelfile, utreefile, ulabelfile, node_origin_file=None, colormap_file=None):
        """
        Creates a Refinement instance from files containing edges, labels, and additional data.

        :param treefile: Path to the file containing edge data.
        :param labelfile: Path to the file containing label data.
        :param utreefile: Path to the file containing unrefined tree data.
        :param ulabelfile: Path to the file containing unrefined label data.
        :param node_origin_file: Optional path to the file containing node origin data.
        :param colormap_file: Optional path to the file containing colormap data.
        :return: A Refinement instance.
        :raises ValueError: If the file format is incorrect.
        """
        with open(treefile, 'r') as f:
            edges = []
            timestamps = {}
            for edge in f:
                try:
                    a, b, t = edge.split()
                    edges.append((a, b))
                    if t != '-1':
                        timestamps[(a,b)] = t
                except:
                    try:
                        a, b = edge.split()
                        edges.append((a, b))
                        timestamps = None
                    except:
                        raise ValueError('Ill-formatted input tree file')
        with open(labelfile, 'r') as f:
            labels = {}
            for line in f:
                a, b = line.split()
                labels[a] = b
        if node_origin_file is not None:
            with open(node_origin_file, 'r') as f:
                node_of_origin = {}
                for line in f:
                    a, b = line.split()
                    node_of_origin[a] = b
        else:
            node_of_origin = None
        return Refinement( edges, labels, 
                          MultiLabeledTree.from_files(utreefile, ulabelfile, colormap_file), 
                          node_of_origin, timestamps=timestamps)
    
    def get_label(self, node):
        """
        Returns the label of a given node.

        :param node: The node whose label is to be retrieved.
        :return: The label of the node.
        """
        return self._tree.nodes[node]['label']

    def _get_comigrations(self):
        """
        Calculates the comigrations and timestamps for the tree.
        """
        import gurobipy as gp
        from gurobipy import GRB
        migs = self.migrations
        migs_st = set([(self.get_label(u), self.get_label(v)) for (u,v) in self.migrations])
        self.X = set()
        for leaf in self.paths:
            path = self.paths[leaf]
            current = None
            for uv in path:
                if uv in migs:
                    if current is None:
                        current = uv
                    else:
                        self.X.add((current, uv))
                        current = uv

        m = gp.Model('PCC')
        m.setParam(GRB.param.LogToConsole, 0)
        m.setParam(GRB.param.LogFile, '')
        m.setParam(GRB.Param.Threads, 0)

        l = m.addVars(migs, range(len(migs)), vtype=GRB.BINARY, lb=0, ub=1)
        pi = m.addVars(range(len(migs)), migs_st, vtype=GRB.CONTINUOUS, lb=0, ub=1)

        for u, v in migs:
            m.addConstr(l.sum(u, v, '*') == 1)
        if len(migs) > 60:
            for E in range(len(migs)):
                for (u, v), (up, vp) in self.X:
                    m.addConstr(l.sum(u, v, range(E)) >= l.sum(up, vp, range(E)))
        else:
            for (u, v), (up, vp) in self.X:
                sum1 = 0
                for E in range(len(migs)):
                    sum1 += 2**(len(migs) - E) * (l[u, v, E] - l[up, vp, E])
                m.addConstr( sum1 >= 0)

        for e in range(len(migs)):
            m.addConstr( pi.sum(e, '*', '*') <= 1 )

        for u, v in migs:
            for e in range(len(migs)):
                m.addConstr(pi[e, self.get_label(u), self.get_label(v)] >= l[u, v, e])

        for e in range(len(migs) - 1):
            m.addConstr(pi.sum(e, '*', '*') >= pi.sum(e + 1, '*', '*'))

        m.setObjective(pi.sum(), GRB.MINIMIZE)
        m.optimize()

        self.timestamps = {}
        self.comigrations = defaultdict(list)
        for u, v in migs:
            for e in range(len(migs)):
                if l[u, v, e].X > 0.5:
                    self.timestamps[(u,v)] = e
                    self.comigrations[e].append((u,v))

    def _greedy_comigrations(self):
        """
        Greedily calculates the comigrations and timestamps for the tree.
        """
        label2comig = defaultdict(list)
        self.comigrations = dict()
        for u, v in self.migrations:
            b = True
            for c in label2comig[(self.get_label(u), self.get_label(v))]:
                bb = True
                for _, vv in self.comigrations[c]:
                    if nx.has_path(self._tree, vv, u):
                        bb = False
                        break
                if bb:
                    self.comigrations[c].append((u,v))
                    self.timestamps[(u,v)] = c
                    b = False
                    break
            if b:
                label2comig[(self.get_label(u), self.get_label(v))].append(len(self.comigrations))
                self.timestamps[(u,v)] = len(self.comigrations)
                self.comigrations[len(self.comigrations)] = [(u,v)]


    def migration_graph(self):
        """
        Returns the migration graph associated with the tree.

        :return: A MigrationGraph object.
        """
        if hasattr(self, '_graph'):
            return self._graph
        self._graph = MigrationGraph(self)
        return self._graph
    
    def migrating_clones(self):
        """
        Returns the set of migrating clones.

        :return: A set of migrating clone nodes.
        """
        return set([self.node_of_origin[u] for u, _ in self.migrations])

    def migration_pattern(self):
        """
        Returns the migration pattern of the tree.

        :return: A tuple (clonality, pattern) where:
                 - clonality is 'm' (monoclonal) or 'p' (polyclonal).
                 - pattern is 'S' (single-source), 'M' (multi-source), or 'R' (reticulate).
        """
        return self.migration_graph().migration_pattern()
    
    def phyleticity(self):
        """
        Determines the phyleticity of the tree.

        :return: 'm' if monoclonal, 'p' if polyclonal.
        """
        mcs = self.migrating_clones()
        for u in mcs:
            is_root_mc = True
            for v in mcs:
                if self.get_path(u, v) is None:
                    is_root_mc = False
            if is_root_mc:
                return 'm'
        return 'p'

    def write_labeling(self, filename):
        """
        Writes the node labels to a file.

        :param filename: Path to the output file.
        """
        with open(filename, 'w') as f:
            for node in self.nodes:
                f.write(f'{node}\t{self.get_label(node)}\n')

    def write_tree_with_timestamps(self, filename):
        """
        Writes the tree edges with timestamps to a file.

        :param filename: Path to the output file.
        """
        with open(filename, 'w+') as f:
            for edge in self.edges:
                if edge in self.migrations:
                    f.write(f'{edge[0]}\t{edge[1]}\t{self.timestamps[edge]}\n')
                else:
                    f.write(f'{edge[0]}\t{edge[1]}\t-1\n')

    def write_node_of_origin(self, filename):
        """
        Writes the node of origin mapping to a file.

        :param filename: Path to the output file.
        """
        with open(filename, 'w') as f:
            for node in self.nodes:
                f.write(f'{node}\t{self.node_of_origin[node]}\n')

    def write_refinement(self, filename_prefix):
        """
        Writes the refined tree, labels, and node origins to files.

        :param filename_prefix: Prefix for the output files.
        """
        self.write_tree_with_timestamps(filename_prefix+'.refined.tree')
        self.write_labeling(filename_prefix+'.location.labeling')
        self.write_node_of_origin(filename_prefix+'.nodes')

    def draw(self, filename=None):
        """
        Draws the Refinement tree using Graphviz.

        :param filename: Optional path to save the rendered graph.
        :return: A Graphviz Digraph object.
        """
        colormap = self.unrefined_tree._colormap
        t = gv.Digraph(node_attr={'penwidth': '3', 'colorscheme': 'set19'}, edge_attr={
                       'penwidth': '3', 'colorscheme': 'set19'})
        t_leaves = gv.Digraph(graph_attr={'rank': 'same'}, node_attr={
                              'penwidth': '3', 'colorscheme': 'set19'})
        for l in self.leaves:
            t_leaves.node(l, f"{l}\n{self.get_label(l)}", color=str(colormap[self.get_label(l)]))
        t.subgraph(t_leaves)
        for v in self.nodes:
            if v not in self.leaves:
                t.node(v, color=str( colormap[self.get_label(v)]))
        for i, j in self.edges:
            if self.get_label(i) != self.get_label(j):
                t.edge( i, j, color=f"{colormap[self.get_label(i)]};0.5:{colormap[self.get_label(j)]}", label=f'  {self.timestamps[(i,j)]}')
            else:
                t.edge( i, j, color=f"{colormap[self.get_label(i)]};0.5:{colormap[self.get_label(j)]}")
        if filename is not None:
            t.render(filename)
        return t


    def n_unobserved_clones(self):
        """
        Returns the number of unobserved clones.

        :return: The number of unobserved clones as an integer.
        """
        return len(self.unobserved_clones)

    def n_migrations(self):
        """
        Returns the number of migrations.

        :return: The number of migrations as an integer.
        """
        return len(self.migrations)

    def n_comigrations(self):
        """
        Returns the number of comigrations.

        :return: The number of comigrations as an integer.
        """
        return len(self.comigrations)

    def n_seeding_locations(self):
        """
        Returns the number of seeding locations.

        :return: The number of seeding locations as an integer.
        """
        return len(self.seeding_locations)

    def parsimony_scores(self):
        """
        Returns the parsimony scores for the tree.

        :return: A dictionary containing the parsimony scores for:
                 - 'U': Unobserved clones.
                 - 'M': Migrations.
                 - 'C': Comigrations.
                 - 'S': Seeding locations.
        """
        return {'U': self.n_unobserved_clones(), 'M': self.n_migrations(),\
                'C': self.n_comigrations(), 'S': self.n_seeding_locations()}