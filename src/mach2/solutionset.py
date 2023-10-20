import networkx as nx
import pandas as pd
from collections import defaultdict
from . import utils


class SolutionSet:

    class Solution:

        def __init__(self, tree, graph):
            self.phylogeny = tree
            self.migration_graph = graph
            self.locations = self.phylogeny.locations
            self.unrefined_phylogeny = tree.unrefined_phylogeny
            self.n_migrations = self.migration_graph.n_mig
            self.n_comigrations = self.migration_graph.n_comig
            self.colormap = self.phylogeny.colormap

    def __init__(self, sol_list):
        self.solution_set = [SolutionSet.Solution(tree, graph) for (tree, graph) in sol_list]
        # self.solution_set = defaultdict(set)
        # for solution in self.solution_set_whole:
        #     self.solution_set[solution.phylogeny.primary_site].add(solution)
        # self.locations = sol_list[0][0].locations

    def __len__(self):
        return len(self.solution_set)
    
    def __getitem__(self, key):
        return self.solution_set[key]
        # if type(key) == type('a'):
        #     return self.solution_set[key]
        # elif type(key) == type(0):
        #     return self.solution_set[key]
        # else:
        #     raise Exception('key is neither an index nor an anatomical site')

    def co_occurence_table(self):
        # if primary is None:
        #     sol_set = self.solution_set_whole
        # else:
        #     sol_set = self.solution_set[primary]
        table = defaultdict(lambda: defaultdict(int))
        for solution in self.solution_set:
            for u1, v1 in solution.migration_graph.migration_edges():
                for u2, v2 in solution.migration_graph.migration_edges():
                    table[f"{u1} -> {v1}"][f"{u2} -> {v2}"] += 1
        return pd.DataFrame(table)
    
    def compute_summary(self, primary=None):
        summary_graph = defaultdict(int)
        for solution in self.solution_set:
            for u1, v1 in solution.migration_graph.migration_edges():
                summary_graph[(u1, v1)] += 1
        return summary_graph

    def summary(self, primary=None, colormap=None, colormap_file=None, consider_multi_edges=False, dot=True):
        # if primary is None:
        #     sol_set = self.solution_set
        # else:
        #     sol_set = self.solution_set
        summary_graph = defaultdict(int)
        locations = set()
        max_val = 0
        for solution in self.solution_set:
            for u1, v1 in solution.migration_graph.migration_edges():
                summary_graph[(u1, v1)] += 1
                if max_val < summary_graph[(u1, v1)]:
                    max_val = summary_graph[(u1, v1)]
                locations.update([u1, v1])

        if dot:
            import graphviz as gv
            if colormap is None:
                if colormap_file is None:
                    if self[0].colormap is None:
                        colormap = utils.get_colormap(self.locations)
                    else:
                        colormap = self[0].colormap
                else:
                    colormap = utils.process_colormap_file(colormap_file)
            dot = gv.Digraph(engine='dot')
            for i in locations:
                dot.node(i, colorscheme='set19', color=str(
                    colormap[i]), shape='box', penwidth='3')

            for uv in summary_graph:
                dot.edge(uv[0], uv[1], penwidth=str(summary_graph[uv] * 5/max_val),
                        colorscheme='set19', color=f'{colormap[uv[0]]};0.5:{colormap[uv[1]]}', label=str(summary_graph[uv]))
            return dot
        else:
            return summary_graph
    
    def summary_dot(self, filename, primary=None, colormap=None, colormap_file=None, consider_multi_edges=False):
        # if primary is None:
        #     sol_set = self.solution_set_whole
        # else:
        #     sol_set = self.solution_set[primary]
        summary_graph = defaultdict(int)
        locations = set()
        max_val = 0
        for solution in self.solution_set:
            for u1, v1 in solution.migration_graph.migration_edges():
                summary_graph[(u1, v1)] += 1
                if max_val < summary_graph[(u1, v1)]:
                    max_val = summary_graph[(u1, v1)]
                locations.update([u1, v1])

        import graphviz as gv
        if colormap is None:
            if colormap_file is None:
                colormap = utils.get_colormap(self.locations)
            else:
                colormap = utils.process_colormap_file(colormap_file)
        node_index = {j:i for i,j in enumerate(locations)}
        with open(filename, 'w+') as f:
            f.write('digraph G {\n')
            for s in locations:
                f.write(f'\t{node_index[s]} [shape=box,penwidth=3,colorscheme=set19,color={colormap[s]},' + f'label="{s}"]\n')
            for st in summary_graph:
                f.write(f'\t{node_index[st[0]]} -> {node_index[st[1]]} [penwidth={summary_graph[st] * 5/max_val},colorscheme=set19,' +
                            f'color="{colormap[st[0]]};0.5:{colormap[st[1]]}"]\n')
            f.write('}\n')

    def write_json(self, name='', filename=None):
        if name == '' and filename is not None:
            name = filename.split('.')[0]
        json_dict = {
                        'name' : name,
                        'original' : {
                            'tree': [[u,v] for u,v in self[0].phylogeny.unrefined_phylogeny.tree.edges], 
                            'labeling': [ [v, self[0].phylogeny.unrefined_phylogeny.get_label(v)] for v in self[0].phylogeny.unrefined_phylogeny.leaves]
                        },
                        'solutions' : [
                            {
                                'name' : f'T-{i}',
                                'tree': [[u,v, (solution.phylogeny.timestamps[(u,v)] if (u,v) in solution.phylogeny.timestamps else -1)] for u,v in solution.phylogeny.tree.edges], 
                                'labeling': [ [v, l['label']] for v, l in solution.phylogeny.tree.nodes.items()], 
                                'migration': [[u,v, solution.migration_graph.n_migrations(u,v)] for (u,v) in solution.migration_graph.migration_edges()],
                                'origin_node': [[ i,v ] for i,v in solution.phylogeny.node_of_origin.items()]
                            } 
                            for i, solution in enumerate(self)
                        ],
                        'summary': [ [u, v, k]  for (u, v), k in self.compute_summary().items()]
                    }
        import json
        if filename is None:
            return json.dumps(json_dict, indent=4)
        else:
            with open(filename, 'w+') as out:
                json.dump(json_dict, out)

    # def summary_write(self, filename, primary=None):
    #     if primary is None:
    #         sol_set = self.solution_set_whole
    #     else:
    #         sol_set = self.solution_set[primary]
    #     summary_graph = defaultdict(int)
    #     locations = set()
    #     max_val = 0
    #     for solution in sol_set:
    #         for u1, v1 in solution.migration_graph.migration_edges():
    #             summary_graph[(u1, v1)] += 1
    #             if max_val < summary_graph[(u1, v1)]:
    #                 max_val = summary_graph[(u1, v1)]
    #             locations.update([u1, v1])

    #     with open(filename, 'w+') as f:
    #         for st in summary_graph:
    #             f.write(f'{st[0]}\t{st[1]}\t{summary_graph[st] * 5/max_val}\n')