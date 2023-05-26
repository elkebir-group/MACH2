from __future__ import annotations
import networkx as nx
import pandas as pd
from . import utils


class Phylogeny:
    """
    Phylogenetic tree with leaf/vertex labeling.

    Attributes:
        tree (nx.DiGraph): Phylogenetic tree with vertex labels.
        colormap (dict[str, int] | None): Mapping of label colors.
        leaves (list[str]): List of leaves.
        nodes (list[str]): List of nodes.
        edges (list[tuple[str, str]]): List of edges.
        root (str): Root node.
        locations (list(str)): Location labels used in the tree.
        n_nodes (int): Number of nodes
        n_edges (int); Number of edges.
        n_leaves (int); Number of leaves.
        n_locations (int); Number of locations.
        vertex_labeled (bool): True if tree is vertex labeled, false if leaf labeled.
        max_height (int): Maximum height of the tree
        paths (dict[str, list[str]]): mapping of each leaf l with the path from root to l.
    """

    def __init__(self, edges: list[tuple[str, str]], labels: dict[str, str], colormap: dict[str, int] | None = None):
        """Initialize the phylogeny from a list of edges and a mapping from leaves/vertices to labels.

        Args:
            edges (list[tuple(str, str)]): list of edges, where each edge is represented by a tuple of two vertices.
            labels (dict[str, str]): Mapping from vertices to labels.
            colormap (dict[str, int] | None): Mapping of label colors (optional). Colors are represented by indices from Brewer color palette set19.

        Raises:
            ValueError: If the input edges do not represent a tree or if the label of any leaf is missing.
        """
        self.tree: nx.DiGraph = nx.DiGraph()
        for edge in edges:
            self.tree.add_edge(edge[0], edge[1])
        if not nx.is_tree(self.tree):
            raise ValueError('Not a valid tree')
        labels = {i: {'label': j} for (i, j) in labels.items()}
        nx.set_node_attributes(self.tree, labels)
        self.colormap: dict[str, int] = colormap

        self.leaves: list[str] = [i for i in self.tree.nodes if self.tree.out_degree(i) == 0]
        self.root: str = [i for i in self.tree.nodes if self.tree.in_degree(i) == 0][0]
        self.nodes: list[str] = [v for v in self.tree.nodes]
        self.edges: list[tuple[str, str]] = [e for e in self.tree.edges]
        self.locations: list[str] = set(s['label'] for s in labels.values())

        self.n_nodes: int = len(self.nodes)
        self.n_edges: int = len(self.edges)
        self.n_leaves: int = len(self.leaves)
        self.n_locations: int = len(self.locations)

        if 'label' in self.tree.nodes[self.root]:
            self.vertex_labeled: bool = True
        else:
            self.vertex_labeled: bool = False

        for leaf in self.leaves:
            if 'label' not in self.tree.nodes[leaf]:
                raise ValueError(f'Site label of leaf {leaf} is missing')
        self._infer_paths()

    def from_file(phylogeny_filename: str, labeling_filename: str) -> Phylogeny:
        """Create a Phylogeny object from tree and labeling files.

        Args:
            phylogeny_filename (str): File name/path of the tree file.
            labeling_filename (str): File name/path of the labeling file.

        Returns:
            Phylogeny: A Phylogeny object created from the provided files.

        Raises:
            ValueError: If the tree file or labeling file is not properly formatted.
        """
        with open(phylogeny_filename, 'r') as phylogeny_file:
            edges = []
            for edge in phylogeny_file:
                try:
                    a, b = edge.split()
                except:
                    raise ValueError('Ill-formatted input tree file')
                edges.append((a, b))
        with open(labeling_filename, 'r') as labeling_file:
            attrs = dict()
            for label in labeling_file:
                try:
                    a, l = label.split()
                except:
                    raise ValueError('Ill-formatted input labeling file')
                attrs[a] = l
        return Phylogeny(edges, attrs)
    
    def from_pandas(phylogeny_df: pd.DataFrame, labeling_df: pd.DataFrame) -> Phylogeny:
        """Create a Phylogeny object from pandas DataFrames containing tree and leaf/vertex labeling.

        Args:
            phylogeny_df (pd.DataFrame): DataFrame representing the tree.
            labeling_df (pd.DataFrame): DataFrame representing the leaf/vertex labeling.

        Returns:
            Phylogeny: A Phylogeny object created from the provided DataFrames.
        """
        return Phylogeny([ (i[0], i[1]) for i in  phylogeny_df.values.tolist()], { i[0]: i[1] for i in  labeling_df.values.tolist()})

    def get_label(self, node: str) -> str:
        """Get the label associated with a given node.

        Args:
            node (str): Node identifier.

        Returns:
            str: Label associated with the node.

        Raises:
            ValueError: If the node does not have a label.
        """
        try:
            return self.tree.nodes[node]['label']
        except:
            raise ValueError(f'Node {node} does not have a label')

    def get_parent_arc(self, node: str) -> tuple[str, str] | None:
        """Get the parent arc of a given node.

        Args:
            node (str): Node identifier.

        Returns:
            tuple[str, str] | None: Parent arc represented as a tuple (parent, node).
                Returns None if the node has no parent, i.e. is the root.

        """
        try:
            return [k for k in self.tree.in_edges(nbunch=node)][0]
        except:
            return None

    def get_children_arcs(self, node: str) -> list[tuple[str, str]]:
        """Get the arcs representing the children of a given node.

        Args:
            node (str): Node identifier.

        Returns:
            list[tuple[str, str]]: list of arcs representing the children of the node.
        """
        return [k for k in self.tree.out_edges(nbunch=node)]

    def _infer_paths(self) -> None:
        """Infer the paths from the root to each leaf in the tree."""

        self.paths: dict[str, list[str]] = {}
        self.max_height: int = 0
        for u in self.leaves:
            self.paths[u] = []
            path_u = nx.shortest_path(self.tree, self.root, u)
            for ii in range(len(path_u) - 1):
                self.paths[u].append((path_u[ii], path_u[ii+1]))
            if len(self.paths[u]) > self.max_height:
                self.max_height = len(self.paths[u])

    def write_tree(self, filename: str) -> None:
        """Write the tree to a file.

        Args:
            filename (str): Name of the file to write the tree.

        Returns:
            None
        """
        with open(filename, 'w+') as f:
            for edge in self.tree.edges:
                f.write(f'{edge[0]} {edge[1]}\n')

    def write_labeling(self, filename: str) -> None:
        """Write the leaf/vertex labeling to a file.

        Args:
            filename (str): Name of the file to write the leaf/vertex labeing

        Returns:
            None
        """
        with open(filename, 'w+') as f:
            if self.vertex_labeled:
                for node in self.nodes:
                    f.write(f'{node} {self.get_label(node)}\n')
            else:
                for leaf in self.leaves:
                    f.write(f'{leaf} {self.get_label(leaf)}\n')

    def write_dot(self, filename: str, colormap: dict[str, int] | None = None, colormap_file: str | None = None) -> None:
        """Write the phylogeny with labels to a DOT file. Colormap priority: input colormap > input colormap file > self.colormap provided in the initializer. If none provided, a random one will be generated.

        Args:
            filename (str): Name of the file to write the DOT representation of the tree.
            colormap (dict[str, int] | None): Mapping of label colors (optional). Colors are represented by indices from Brewer color palette set19.
            colormap_file (str | None): File name/path of a colormap file.

        Returns:
            None
        """
        if colormap is None:
            if colormap_file is None:
                if self.colormap is None:
                    colormap = utils.get_colormap(self.locations)
                else:
                    colormap = self.colormap
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

    def draw(self, colormap: dict[str, int] | None = None, colormap_file: str | None = None) -> None:
        """Returns tree as a Graphviz DiGraph. Can be used to draw in JupyterLab. Colormap priority: input colormap > input colormap file > self.colormap provided in the initializer. If none provided, a random one will be generated.

        Args:
            filename (str): Name of the file to write the DOT representation of the tree.
            colormap (dict[str, int] | None): Mapping of label colors (optional). Colors are represented by indices from Brewer color palette set19.
            colormap_file (str | None): File name/path of a colormap file.

        Returns:
            None
        """
        import graphviz as gv
        if colormap is None:
            if colormap_file is None:
                colormap = utils.get_colormap(self.locations)
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
                    t.node(v, color=str(len(self.locations) + 1))
        for i, j in self.tree.edges:
            if self.vertex_labeled == True:
                t.edge(
                    i, j, color=f"{colormap[self.get_label(i)]};0.5:{colormap[self.get_label(j)]}")
            else:
                t.edge(i, j, color=str(len(self.locations) + 1))
        return t
    
    def _repr_mimebundle_(self, include=None, exclude=None, **_):
        return self.draw()._repr_mimebundle_(include=include, exclude=exclude)

class RefinedPhylogeny(Phylogeny):

    def __init__(self, phylogeny: Phylogeny, raw: dict):
        """Initialize the refined phylogeny from raw output returned by MACH.solve().

        Args:
            phylogeny (Phylogeny): Phylogeny before refinement
            raw (dict): Raw solutions returned by MACH.solve()
        """
        self.unrefined_phylogeny: Phylogeny = phylogeny
        self.n_migrations: int = raw['n_mig']
        self.n_comigrations: int = raw['n_comig']
        ell = raw['vertex_multilabeling']
        barG = raw['aug_mig_graph']
        edges = []
        labels = {}
        self.node_of_origin: dict[str, str] = {}
        for v in phylogeny.nodes:
            if len(ell[v]) > 1:
                for s_1, s_2 in barG[v]:
                    labels[f'{v}^{s_1}'] = s_1
                    labels[f'{v}^{s_2}'] = s_2
                    self.node_of_origin[f'{v}^{s_1}'] = v
                    self.node_of_origin[f'{v}^{s_2}'] = v
                    edges.append((f'{v}^{s_1}', f'{v}^{s_2}'))
            else:
                labels[v] = ell[v][0]
                self.node_of_origin[v] = v
        for u, v in phylogeny.edges:
            s_1, s_2 = barG[(u, v)][0]
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
        self.primary_site: str = self.get_label(self.root)
        self.infer_timestamps()

    def _infer_timestamps_rec(self, u, prev_st, comig_to_mig, mig_to_comig, comig_graph, last_comig):
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
                self._infer_timestamps_rec(v, prev_st, comig_to_mig, mig_to_comig, comig_graph, new_comig)
                prev_st.pop()
            else:
                self._infer_timestamps_rec(v, prev_st, comig_to_mig, mig_to_comig, comig_graph, last_comig)

    def infer_timestamps(self):
        prev_st = []
        comig_to_mig = {}
        mig_to_comig = {}
        comig_graph = nx.DiGraph()
        self._infer_timestamps_rec(self.root, prev_st, comig_to_mig, mig_to_comig, comig_graph, None)
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
                colormap = utils.get_colormap(self.locations)
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
                colormap = utils.get_colormap(self.locations)
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