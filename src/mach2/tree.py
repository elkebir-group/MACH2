import networkx as nx
import pandas as pd
from .migrationgraph import MigrationGraph
from . import utils
import gurobipy as gp
from gurobipy import GRB
from collections import defaultdict


class Tree:
    def __init__(self, edges):
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
        return set(self.edges) == set(other.edges)
    
    def __ne__(self, other):
        return set(self.edges) != set(other.edges)
    
    def __hash__(self):
        return hash(frozenset(self.edges))

    def from_file(filename):
        with open(filename, 'r') as f:
            edges = []
            for edge in f.readlines():
                try:
                    a, b = edge.split()
                except:
                    raise ValueError('Ill-formatted input tree file')
                edges.append((a, b))
        return Tree(edges)
    
    def from_pandas(dataframe):
        pass
    
    def get_parent(self, node):
        try:
            return [k for k in self._tree.in_edges(nbunch=node)][0][0]
        except:
            return None
    
    def get_children(self, node):
        return [k[1] for k in self._tree.out_edges(nbunch=node)]
    
    def get_descendant_leaves(self, node):
        if node in self.leaves:
            return [node]
        desc = []
        for v in self.get_children(node):
            desc += self.get_descendant_leaves(v)
        return desc
    
    def get_descendants(self, node):
        desc = []
        for v in self.get_children(node):
            desc += [v] + self.get_descendants(v)
        return desc
    
    def get_path(self, node1, node2):
        try:
            return nx.shortest_path(self._tree, node1, node2)
        except:
            return None

    def write_tree(self, filename):
        with open(filename, 'w') as f:
            for edge in self.edges:
                f.write(f'{edge[0]}\t{edge[1]}\n')

    def draw(self):
        import graphviz as gv
        t = gv.Digraph(node_attr={'penwidth': '3'}, edge_attr={'penwidth': '3'})
        t_leaves = gv.Digraph(graph_attr={'rank': 'same'}, node_attr={'penwidth': '3'})
        for l in self.leaves:
            t_leaves.node(l, f"{l}")
        t.subgraph(t_leaves)
        for i, j in self.edges:
            t.edge(i, j)
        return t
    
    def _repr_mimebundle_(self, include=None, exclude=None, **_):
        return self.draw()._repr_mimebundle_(include=include, exclude=exclude)
    

class MultiLabeledTree(Tree):
    def __init__(self, edges, labels, colormap=None):
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
        return super().__eq__(other) and set([(u, frozenset(self.get_labels(u))) for u in self.nodes]) == set([(u, frozenset(other.get_labels(u))) for u in other.nodes])
    
    def __ne__(self, other):
        return not self == other
    
    def __hash__(self):
        return hash((frozenset(self.edges), frozenset([(u, frozenset(self.get_labels(u))) for u in self.nodes])))

    def from_files(treefile, labelfile, colormap_file=None):
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
        return self._tree.nodes[node]['label']
    
    def add_label(self, node, label):
        if label not in self.locations:
            self.locations.add(label)
        self._tree.nodes[node]['label'].add(label)

    def write_labelings(self, filename):
        with open(filename, 'w') as f:
            for node in self.nodes:
                labels = self.get_labels(node)
                if len(labels) > 0:
                    f.write(f'{node}')
                    for label in labels:
                        f.write(f'\t{label}')
                    f.write('\n')

    def set_colormap(self, colormap):
        self._colormap = colormap

    def set_time_of_observation(self, time_of_observation):
        self.time_of_observation = time_of_observation

    def draw(self):
        import graphviz as gv
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
    def __init__(self, edges, labels, unrefined_tree, node_of_origin=None, timestamps=None):
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
            self._get_comigrations()
        else:
            self.timestamps = timestamps
            self.comigrations = defaultdict(list)
            for uv, t in timestamps.items():
                self.comigrations[t].append(uv)
        self.seeding_locations = list(set([attrs[u]['label'] for (u, _) in self.migrations]))

    def __eq__(self, other):
        return super().__eq__(other) and set([(u, self.get_label(u)) for u in self.nodes]) == set([(u, other.get_label(u)) for u in other.nodes])
    
    def __ne__(self, other):
        return not self == other
    
    def __hash__(self):
        return hash((frozenset(self.edges), frozenset([(u, self.get_label(u)) for u in self.nodes])))

    def from_refinement_graph(unrefined_tree, refinement_graph):
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
        return self._tree.nodes[node]['label']

    def _get_comigrations(self):
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

    def migration_graph(self):
        if hasattr(self, '_graph'):
            return self._graph
        self._graph = MigrationGraph(self)
        return self._graph
    
    def migrating_clones(self):
        return set([self.node_of_origin[u] for u, _ in self.migrations])

    def migration_pattern(self):
        return self.migration_graph().migration_pattern()
    
    def phyleticity(self):
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
        with open(filename, 'w') as f:
            for node in self.nodes:
                f.write(f'{node}\t{self.get_label(node)}\n')

    def write_tree_with_timestamps(self, filename):
        with open(filename, 'w+') as f:
            for edge in self.edges:
                if edge in self.migrations:
                    f.write(f'{edge[0]}\t{edge[1]}\t{self.timestamps[edge]}\n')
                else:
                    f.write(f'{edge[0]}\t{edge[1]}\t-1\n')

    def write_node_of_origin(self, filename):
        with open(filename, 'w') as f:
            for node in self.nodes:
                f.write(f'{node}\t{self.node_of_origin[node]}\n')

    def write_refinement(self, filename_prefix):
        self.write_tree_with_timestamps(filename_prefix+'.refined.tree')
        self.write_labeling(filename_prefix+'.location.labeling')
        self.write_node_of_origin(filename_prefix+'.nodes')

    def draw(self, filename=None):
        import graphviz as gv
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
        return len(self.unobserved_clones)

    def n_migrations(self):
        return len(self.migrations)

    def n_comigrations(self):
        return len(self.comigrations)

    def n_seeding_locations(self):
        return len(self.seeding_locations)

    def parsimony_scores(self):
        return {'U': self.n_unobserved_clones(), 'M': self.n_migrations(),\
                'C': self.n_comigrations(), 'S': self.n_seeding_locations()}