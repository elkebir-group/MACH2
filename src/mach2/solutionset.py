import networkx as nx
import pandas as pd
from collections import defaultdict
from . import utils


class SolutionSet:

    class Solution:

        def __init__(self, tree, graph):
            self.clone_tree = tree
            self.migration_graph = graph
            self.sites = self.clone_tree.sites
            self.n_migrations = self.migration_graph.n_mig
            self.n_comigrations = self.migration_graph.n_comig

    def __init__(self, sol_list):
        self.solution_set_whole = [SolutionSet.Solution(tree, graph) for (tree, graph) in sol_list]
        self.solution_set = defaultdict(set)
        for solution in self.solution_set_whole:
            self.solution_set[solution.clone_tree.primary_site].add(solution)
        self.sites = sol_list[0][0].sites

    def __len__(self):
        return len(self.solution_set_whole)
    
    def __getitem__(self, key):
        if type(key) == type('a'):
            return self.solution_set[key]
        elif type(key) == type(0):
            return self.solution_set_whole[key]
        else:
            raise Exception('key is neither an index nor an anatomical site')

    def co_occurence_table(self, primary=None):
        if primary is None:
            sol_set = self.solution_set_whole
        else:
            sol_set = self.solution_set[primary]
        table = defaultdict(lambda: defaultdict(int))
        for solution in sol_set:
            for u1, v1 in solution.migration_graph.migration_edges():
                for u2, v2 in solution.migration_graph.migration_edges():
                    table[f"{u1} -> {v1}"][f"{u2} -> {v2}"] += 1
        return pd.DataFrame(table)

    def summary(self, primary=None, colormap=None, colormap_file=None, consider_multi_edges=False):
        if primary is None:
            sol_set = self.solution_set_whole
        else:
            sol_set = self.solution_set[primary]
        summary_graph = defaultdict(int)
        sites = set()
        max_val = 0
        for solution in sol_set:
            for u1, v1 in solution.migration_graph.migration_edges():
                summary_graph[(u1, v1)] += 1
                if max_val < summary_graph[(u1, v1)]:
                    max_val = summary_graph[(u1, v1)]
                sites.update([u1, v1])

        import graphviz as gv
        if colormap is None:
            if colormap_file is None:
                colormap = utils.get_colormap(self.sites)
            else:
                colormap = utils.process_colormap_file(colormap_file)
        dot = gv.Digraph(engine='dot')
        for i in sites:
            dot.node(i, colorscheme='set19', color=str(
                colormap[i]), shape='box', penwidth='3')

        for uv in summary_graph:
            dot.edge(uv[0], uv[1], penwidth=str(summary_graph[uv] * 5/max_val),
                     colorscheme='set19', color=f'{colormap[uv[0]]};0.5:{colormap[uv[1]]}')
            
        return dot
    
    def summary_dot(self, filename, primary=None, colormap=None, colormap_file=None, consider_multi_edges=False):
        if primary is None:
            sol_set = self.solution_set_whole
        else:
            sol_set = self.solution_set[primary]
        summary_graph = defaultdict(int)
        sites = set()
        max_val = 0
        for solution in sol_set:
            for u1, v1 in solution.migration_graph.migration_edges():
                summary_graph[(u1, v1)] += 1
                if max_val < summary_graph[(u1, v1)]:
                    max_val = summary_graph[(u1, v1)]
                sites.update([u1, v1])

        import graphviz as gv
        if colormap is None:
            if colormap_file is None:
                colormap = utils.get_colormap(self.sites)
            else:
                colormap = utils.process_colormap_file(colormap_file)
        node_index = {j:i for i,j in enumerate(sites)}
        with open(filename, 'w+') as f:
            f.write('digraph G {\n')
            for s in sites:
                f.write(f'\t{node_index[s]} [shape=box,penwidth=3,colorscheme=set19,color={colormap[s]},' + f'label="{s}"]\n')
            for st in summary_graph:
                f.write(f'\t{node_index[st[0]]} -> {node_index[st[1]]} [penwidth={summary_graph[st] * 5/max_val},colorscheme=set19,' +
                            f'color="{colormap[st[0]]};0.5:{colormap[st[1]]}"]\n')
            f.write('}\n')

    # def summary_write(self, filename, primary=None):
    #     if primary is None:
    #         sol_set = self.solution_set_whole
    #     else:
    #         sol_set = self.solution_set[primary]
    #     summary_graph = defaultdict(int)
    #     sites = set()
    #     max_val = 0
    #     for solution in sol_set:
    #         for u1, v1 in solution.migration_graph.migration_edges():
    #             summary_graph[(u1, v1)] += 1
    #             if max_val < summary_graph[(u1, v1)]:
    #                 max_val = summary_graph[(u1, v1)]
    #             sites.update([u1, v1])

    #     with open(filename, 'w+') as f:
    #         for st in summary_graph:
    #             f.write(f'{st[0]}\t{st[1]}\t{summary_graph[st] * 5/max_val}\n')