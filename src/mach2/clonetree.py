import networkx as nx
from . import utils


class CloneTree:

    def __init__(self, edges, labels, infer_paths=False):
        self.tree = nx.DiGraph()
        for edge in edges:
            self.tree.add_edge(edge[0], edge[1])
        if not nx.is_tree(self.tree):
            raise RuntimeError('Not a valid tree')
        labels = {i: {'label': j} for (i, j) in labels.items()}
        nx.set_node_attributes(self.tree, labels)

        self.leaves = [i for i in self.tree.nodes if self.tree.out_degree(i) == 0]
        self.root = [i for i in self.tree.nodes if self.tree.in_degree(i) == 0][0]
        self.nodes = [v for v in self.tree.nodes]
        self.edges = [e for e in self.tree.edges]
        self.sites = set(s['label'] for s in labels.values())

        self.n_nodes = len(self.nodes)
        self.n_edges = len(self.edges)
        self.n_leaves = len(self.leaves)
        self.n_sites = len(self.sites)

        if 'label' in self.tree.nodes[self.root]:
            self.vertex_labeled = True
        else:
            self.vertex_labeled = False

        for leaf in self.leaves:
            if 'label' not in self.tree.nodes[leaf]:
                raise ValueError(f'Site label of clone {leaf} is missing')
        if infer_paths:
            self.infer_paths()

    def from_file(clone_tree_filename, leaf_labeling_filename, infer_paths=False):
        with open(clone_tree_filename, 'r') as clone_tree_file:
            edges = []
            for edge in clone_tree_file:
                try:
                    a, b = edge.split()
                except:
                    raise ValueError('Ill-formatted input clone tree file')
                edges.append((a, b))
        with open(leaf_labeling_filename, 'r') as leaf_labeling_file:
            attrs = dict()
            for label in leaf_labeling_file:
                try:
                    a, l = label.split()
                except:
                    raise ValueError('Ill-formatted input leaf labeling file')
                attrs[a] = l
        return CloneTree(edges, attrs, infer_paths=infer_paths)

    def get_label(self, node):
        try:
            return self.tree.nodes[node]['label']
        except:
            raise ValueError(f'Node {node} does not have a label')

    def get_parent_arc(self, node):
        try:
            return [k for k in self.tree.in_edges(nbunch=node)][0]
        except:
            return None

    def get_children_arcs(self, node):
        return [k for k in self.tree.out_edges(nbunch=node)]

    def infer_paths(self):
        self.paths = {}
        self.max_height = 0
        for u in self.leaves:
            self.paths[u] = []
            path_u = nx.shortest_path(self.tree, self.root, u)
            for ii in range(len(path_u) - 1):
                self.paths[u].append((path_u[ii], path_u[ii+1]))
            if len(self.paths[u]) > self.max_height:
                self.max_height = len(self.paths[u])

    def write_tree(self, filename):
        with open(filename, 'w+') as f:
            for edge in self.tree.edges:
                f.write(f'{edge[0]} {edge[1]}\n')

    def write_labeling(self, filename):
        with open(filename, 'w+') as f:
            if self.vertex_labeled:
                for node in self.nodes:
                    f.write(f'{node} {self.get_label(node)}\n')
            else:
                for leaf in self.leaves:
                    f.write(f'{leaf} {self.get_label(leaf)}\n')

    def write_dot(self, filename, colormap=None, colormap_file=None):
        if colormap is None:
            if colormap_file is None:
                colormap = utils.get_colormap(self.sites)
            else:
                colormap = utils.process_colormap_file(colormap_file)
        with open(filename, 'w+') as f:
            f.write('digraph T {\n\t{\n\t\trank=same\n')
            node_edge_index = {}
            ind = 0
            for l in self.leaves:
                f.write(
                    f'\t\t{ind} [penwidth=3,colorscheme=set19,color={colormap[self.get_label(l)]},label="{l}\\n{self.get_label(l)}"]\n')
                node_edge_index[l] = ind
                ind += 1
            f.write('\t}\n')
            for v in self.tree.nodes:
                if v not in self.leaves:
                    if self.vertex_labeled == True:
                        f.write(
                            f'\t{ind} [penwidth=3,colorscheme=set19,color={colormap[self.get_label(v)]},label="{v}"]\n')
                    else:
                        f.write(
                            f'\t{ind} [penwidth=3,colorscheme=set19,color=0,label="{v}"]\n')
                    node_edge_index[v] = ind
                    ind += 1
            for i, j in self.tree.edges:
                if self.vertex_labeled == True:
                    f.write(f'\t{node_edge_index[i]} -> {node_edge_index[j]} [penwidth=3,colorscheme=set19,' +
                            f'color=\"{colormap[self.get_label(i)]};0.5:{colormap[self.get_label(j)]}\"]\n')
                else:
                    f.write(f'\t{node_edge_index[i]} -> {node_edge_index[j]} [penwidth=3,colorscheme=set19,' +
                            f'color=0]\n')
            f.write('}\n')

    def draw(self, colormap=None, colormap_file=None):
        import graphviz as gv
        if colormap is None:
            if colormap_file is None:
                colormap = utils.get_colormap(self.sites)
            else:
                colormap = utils.process_colormap_file(colormap_file)
        t = gv.Digraph(node_attr={'penwidth': '3', 'colorscheme': 'set19'}, edge_attr={
                       'penwidth': '3', 'colorscheme': 'set19'})
        t_leaves = gv.Digraph(graph_attr={'rank': 'same'}, node_attr={
                              'penwidth': '3', 'colorscheme': 'set19'})
        for l in self.leaves:
            t_leaves.node(l, f"{l}\n{self.get_label(l)}",
                          color=str(colormap[self.get_label(l)]))
        t.subgraph(t_leaves)
        for v in self.tree.nodes:
            if v not in self.leaves:
                if self.vertex_labeled == True:
                    t.node(v, color=str(
                        colormap[self.get_label(v)]))
                else:
                    t.node(v, color='0')
        for i, j in self.tree.edges:
            if self.vertex_labeled == True:
                t.edge(
                    i, j, color=f"{colormap[self.get_label(i)]};0.5:{colormap[self.get_label(j)]}")
            else:
                t.edge(i, j, color='0')
        return t
    
    # def _repr_pretty_(self):
    #     self.draw()

class RefinedCloneTree(CloneTree):

    def __init__(self, clone_tree, raw):
        ell = raw['vertex_multilabeling']
        G = raw['aug_mig_graph']
        edges = []
        labels = {}
        for v in clone_tree.nodes:
            if len(ell[v]) > 1:
                for s_1, s_2 in G[v]:
                    labels[f'{v}^{s_1}'] = s_1
                    labels[f'{v}^{s_2}'] = s_2
                    edges.append((f'{v}^{s_1}', f'{v}^{s_2}'))
            else:
                labels[v] = ell[v][0]
        for u, v in clone_tree.edges:
            s_1, s_2 = G[(u, v)][0]
            if len(ell[u]) > 1:
                edge = (f'{u}^{s_1}',)
            else:
                edge = (u,)
            if len(ell[v]) > 1:
                edge += (f'{v}^{s_2}',)
            else:
                edge += (v,)
            edges.append(edge)
        super().__init__(edges, labels)
        self.primary_site = self.get_label(self.root)
        self.infer_timestamps()

    def infer_timestamps_rec(self, u, prev_st, comig_to_mig, mig_to_comig, comig_graph, last_comig):
        s = self.get_label(u)
        for _, v in self.get_children_arcs(u):
            t = self.get_label(v)
            if s != t:
                if (s,t,1) not in comig_to_mig:
                    comig_to_mig[(s, t, 1)] = [(u,v)]
                    mig_to_comig[(u,v)] = (s, t, 1)
                    if last_comig is not None:
                        comig_graph.add_edge(last_comig, (s, t, 1))
                    else:
                        comig_graph.add_node((s, t, 1))
                    new_comig = (s,t,1)
                elif (s,t) in prev_st:
                    l = len([i for i in prev_st if (i[0]==s and i[1]==t)])
                    if (s,t,l+1) in comig_to_mig:
                        comig_to_mig[(s, t, l+1)].append((u,v))
                        mig_to_comig[(u,v)] = (s, t, l+1)
                        if last_comig is not None:
                            comig_graph.add_edge(last_comig, (s, t, l+1))
                        else:
                            comig_graph.add_node((s, t, l+1))
                        new_comig = (s,t,l+1)
                    else:
                        comig_to_mig[(s, t, l+1)] = [(u,v)]
                        mig_to_comig[(u,v)] = (s, t, 1)
                        if last_comig is not None:
                            comig_graph.add_edge(last_comig, (s, t, 1))
                        else:
                            comig_graph.add_node((s, t, 1))
                        new_comig = (s,t,1)
                else:
                    comig_to_mig[(s, t, 1)].append((u,v))
                    mig_to_comig[(u,v)] = (s, t, 1)
                    if last_comig is not None:
                        comig_graph.add_edge(last_comig, (s, t, 1))
                    else:
                        comig_graph.add_node((s, t, 1))
                    new_comig = (s,t,1)
                prev_st.append((s,t))
                self.infer_timestamps_rec(v, prev_st, comig_to_mig, mig_to_comig, comig_graph, new_comig)
                prev_st.pop()
            else:
                self.infer_timestamps_rec(v, prev_st, comig_to_mig, mig_to_comig, comig_graph, last_comig)

    def infer_timestamps(self):
        prev_st = []
        comig_to_mig = {}
        mig_to_comig = {}
        comig_graph = nx.DiGraph()
        self.infer_timestamps_rec(self.root, prev_st, comig_to_mig, mig_to_comig, comig_graph, None)
        comig_to_ts = {j:i+1 for i,j in enumerate(nx.topological_sort(comig_graph))}
        timestamps = {}
        for u, v in self.edges:
            if self.get_label(u) != self.get_label(v):
                timestamps[(u,v)] = comig_to_ts[mig_to_comig[(u,v)]]
        self.timestamps = timestamps

    def write_tree(self, filename):
        with open(filename, 'w+') as f:
            for edge in self.tree.edges:
                if self.get_label(edge[0]) != self.get_label(edge[1]):
                    f.write(f'{edge[0]} {edge[1]} {self.timestamps[edge]}\n')
                else:
                    f.write(f'{edge[0]} {edge[1]} -1\n')

    def write_labeling(self, filename):
        with open(filename, 'w+') as f:
            for node in self.nodes:
                f.write(f'{node} {self.get_label(node)}\n')

    def write_dot(self, filename, colormap=None, colormap_file=None):
        if colormap is None:
            if colormap_file is None:
                colormap = utils.get_colormap(self.sites)
            else:
                colormap = utils.process_colormap_file(colormap_file)
        with open(filename, 'w+') as f:
            f.write('digraph T {\n\t{\n\t\trank=same\n')
            node_edge_index = {}
            ind = 0
            for l in self.leaves:
                f.write(
                    f'\t\t{ind} [penwidth=3,colorscheme=set19,color={colormap[self.get_label(l)]},label="{l}\\n{self.get_label(l)}"]\n')
                node_edge_index[l] = ind
                ind += 1
            f.write('\t}\n')
            for v in self.tree.nodes:
                if v not in self.leaves:
                    f.write(
                        f'\t{ind} [penwidth=3,colorscheme=set19,color={colormap[self.get_label(v)]},label="{v}"]\n')
                    node_edge_index[v] = ind
                    ind += 1
            for i, j in self.tree.edges:
                if self.get_label(i) != self.get_label(j):
                    f.write(f'\t{node_edge_index[i]} -> {node_edge_index[j]} [penwidth=3,colorscheme=set19,' +
                        f'color=\"{colormap[self.get_label(i)]};0.5:{colormap[self.get_label(j)]}\",label=\"  {self.timestamps[(i,j)]}\"]\n')
                else:
                    f.write(f'\t{node_edge_index[i]} -> {node_edge_index[j]} [penwidth=3,colorscheme=set19,' +
                        f'color=\"{colormap[self.get_label(i)]};0.5:{colormap[self.get_label(j)]}\"]\n')
            f.write('}\n')

    def draw(self, colormap=None, colormap_file=None):
        import graphviz as gv
        if colormap is None:
            if colormap_file is None:
                colormap = utils.get_colormap(self.sites)
            else:
                colormap = utils.process_colormap_file(colormap_file)
        t = gv.Digraph(node_attr={'penwidth': '3', 'colorscheme': 'set19'}, edge_attr={
                       'penwidth': '3', 'colorscheme': 'set19'})
        t_leaves = gv.Digraph(graph_attr={'rank': 'same'}, node_attr={
                              'penwidth': '3', 'colorscheme': 'set19'})
        for l in self.leaves:
            t_leaves.node(l, f"{l}\n{self.get_label(l)}",
                          color=str(colormap[self.get_label(l)]))
        t.subgraph(t_leaves)
        for v in self.tree.nodes:
            if v not in self.leaves:
                t.node(v, color=str(
                    colormap[self.get_label(v)]))
        for i, j in self.tree.edges:
            if self.get_label(i) != self.get_label(j):
                t.edge(
                    i, j, color=f"{colormap[self.get_label(i)]};0.5:{colormap[self.get_label(j)]}", label=f'  {self.timestamps[(i,j)]}')
            else:
                t.edge(
                    i, j, color=f"{colormap[self.get_label(i)]};0.5:{colormap[self.get_label(j)]}")
        return t