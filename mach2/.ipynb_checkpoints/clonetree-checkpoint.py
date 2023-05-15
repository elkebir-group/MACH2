import networkx as nx

class CloneTree:

    def __init__(self, edges, labels, primary_site = None):
        self.tree = nx.DiGraph()
        for edge in edges:
            self.tree.add_edge(edge[0], edge[1])
        if not nx.is_tree(self.tree):
            raise RuntimeError('Not a valid tree')
        labels = { i: {'label': j} for (i, j) in labels.items()}
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

        for leaf in self.leaves:
            if 'label' not in self.tree.nodes[leaf]:
                raise ValueError(f'Site label of clone {leaf} is missing')
            if primary_site is not None and primary_site not in self.sites:
                raise ValueError('Primary site not found')
            else:
                self.primary_site = primary_site

    def from_file(clone_tree_filename, leaf_labeling_filename, primary_site = None):
        with open(clone_tree_filename, 'r') as clone_tree_file:
            edges = []
            for edge in clone_tree_file:
                try:
                    a, b = edge.split()
                except:
                    raise ValueError('Ill-formatted input clone tree file')
                edges.append((a,b))
        with open(leaf_labeling_filename, 'r') as leaf_labeling_file:
            attrs = dict()
            for label in leaf_labeling_file:
                try:
                    a, l = label.split()
                except:
                    raise ValueError('Ill-formatted input leaf labeling file')
                attrs[a] = l
        return CloneTree(edges, attrs, primary_site)
    
    def get_label(self, node):
        return self.tree.nodes[node]['label']

    def get_parent_arc(self, node):
        return [k for k in self.tree.in_edges(nbunch=node)][0]

    def get_children_arcs(self, node):
        return [k for k in self.tree.out_edges(nbunch=node)]

    def infer_paths(self):
        self.paths = {}
        self.max_height = 0
        for u in self.leaves:
            self.paths[u] = []
            path_u = nx.shortest_path(self.tree, self.root, u)
            for ii in range(len(path_u) - 1):
                self.paths[u].append( (path_u[ii], path_u[ii+1]) )
            if len(self.paths[u]) > self.max_height:
                self.max_height = len(self.paths[u])

    def write_dot(self, filename, colormap=None, vertex_labels=True):
        with open(filename, 'w+') as f:
            f.write('digraph T {\n\t{\n\t\trank=same\n')
            node_edge_index = {}
            ind = 0
            for l in self.leaves:
                f.write(f'\t\t{ind} [penwidth=3,colorscheme=set19,color={colormap[self.get_label(l)]},label="{l}\\n{self.get_label(l)}"]\n')
                node_edge_index[l] = ind
                ind += 1
            f.write('\t}\n')
            for v in self.tree.nodes:
                if v not in self.leaves:
                    if vertex_labels == True:
                        f.write(f'\t{ind} [penwidth=3,colorscheme=set19,color={colormap[self.get_label(v)]},label="{v}"]\n')
                    else:
                        f.write(f'\t{ind} [penwidth=3,colorscheme=set19,color=0,label="{v}"]\n')
                    node_edge_index[v] = ind
                    ind += 1
            for i,j in self.tree.edges:
                if vertex_labels == True:
                    f.write(f'\t{node_edge_index[i]} -> {node_edge_index[j]} [penwidth=3,colorscheme=set19,' + 
                        f'color=\"{colormap[self.get_label(i)]};0.5:{colormap[self.get_label(j)]}\"]\n')
                else:
                    f.write(f'\t{node_edge_index[i]} -> {node_edge_index[j]} [penwidth=3,colorscheme=set19,' + 
                        f'color=0]\n')
            f.write('}\n')

    def write_tree(self, filename):
        with open(filename, 'w+') as f:
            for edge in self.tree.edges:
                f.write(f'{edge[0]} {edge[1]}\n')

    def write_labeling(self, filename):
        with open(filename, 'w+') as f:
            for node in self.tree.nodes:
                f.write(f'{node} {self.get_label(node)}\n')

    def draw(self):
        pass

