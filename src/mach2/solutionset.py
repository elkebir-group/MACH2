import networkx as nx
import pandas as pd
from collections import defaultdict
import glob
from . import utils
from .tree import Refinement


class SolutionSet:

    def __init__(self, sol_list=[], check=True):
        self.solution_set = set(sol_list)
        self.unrefined_tree = [i for i in sol_list][0].unrefined_tree if len(sol_list) > 0 else None
        if check:
            for sol in self.solution_set:
                if self.unrefined_tree != sol.unrefined_tree:
                    raise ValueError('Solutions with different original unrefined trees')
                
    def from_files(output_prefix, utreefile, ulabelfile, colormap_file=None):
        sol_list = []
        for treefile in glob.glob(output_prefix + '*.refined.tree'):
            sol_list.append(Refinement.from_files(treefile, 
                                                  treefile[:-13]+'.location.labeling',
                                                  utreefile, ulabelfile,
                                                  treefile[:-13]+'.nodes',
                                                  colormap_file))
        return SolutionSet(sol_list, check=False)

    def __len__(self):
        return len(self.solution_set)
    
    def __iter__(self):
        return iter(self.solution_set)
    
    def __add__(s1, s2):
        if len(s1.solution_set) == 0:
            return s2
        elif len(s2.solution_set) == 0:
            return s1
        return SolutionSet(s1.solution_set.union(s2.solution_set), check=(s1.unrefined_tree!=s2.unrefined_tree))
    
    def filter(self, criteria_ordering):
        if criteria_ordering is not None:
            sort_key_map = { 
                'M': lambda ref: len(ref.migrations), 
                'U': lambda ref: len(ref.unobserved_clones), 
                'C': lambda ref: len(ref.comigrations) 
            }
            sorter = lambda ref, co: tuple(sort_key_map[letter](ref) for letter in co)
            min_val = min([sorter(sol, criteria_ordering) for sol in self.solution_set])
            return SolutionSet([sol for sol in self.solution_set if sorter(sol, criteria_ordering) == min_val], check=False)


    # def co_occurence_table(self):
    #     table = defaultdict(lambda: defaultdict(int))
    #     for solution in self.solution_set:
    #         for u1, v1 in solution.migration_graph.migration_edges():
    #             for u2, v2 in solution.migration_graph.migration_edges():
    #                 table[f"{u1} -> {v1}"][f"{u2} -> {v2}"] += 1
    #     return pd.DataFrame(table)
    
    # def compute_summary(self, primary=None):
    #     summary_graph = defaultdict(int)
    #     for solution in self.solution_set:
    #         for u1, v1 in solution.migration_graph.migration_edges():
    #             summary_graph[(u1, v1)] += 1
    #     return summary_graph

    def summary(self, draw=True):
        summary_graph = defaultdict(int)
        locations = set()
        max_val = 0
        for solution in self.solution_set:
            for u1, v1 in solution.migration_graph().get_migrations():
                summary_graph[(u1, v1)] += 1
                if max_val < summary_graph[(u1, v1)]:
                    max_val = summary_graph[(u1, v1)]
                locations.update([u1, v1])

        if draw:
            import graphviz as gv
            colormap = self.unrefined_tree._colormap
            dot = gv.Digraph(engine='dot')
            for i in locations:
                dot.node(i, colorscheme='set19', color=str(colormap[i]), shape='box', penwidth='3')

            for uv in summary_graph:
                dot.edge(uv[0], uv[1], penwidth=str(7) if summary_graph[uv] == max_val else str(3),
                        colorscheme='set19', color=f'0,0,{1 - summary_graph[uv]/max_val}', 
                        label=f'{summary_graph[uv]}/{len(self)}')
            return dot
        else:
            return summary_graph
        
    def write(self, output_prefix):
        for i, sol in enumerate(self):
            sol.write_refinement(f'{output_prefix}sol{i}')

    def json(self, name=''):
        return {
                    'name' : name,
                    'original' : {
                        'tree': [[u,v] for u,v in self.unrefined_tree.edges], 
                        'labeling': [ [v, list(self.unrefined_tree.get_labels(v))] for v in self.unrefined_tree.nodes]
                    },
                    'solutions' : [
                        {
                            'name' : f'\#{i}',
                            'tree': [[u,v, (sol.timestamps[(u,v)] if (u,v) in sol.timestamps else -1)] for u,v in sol.edges], 
                            'labeling': [ [v, sol.get_label(v)] for v in sol.nodes], 
                            'migration': [[s, t, c] for (s,t), c in sol.migration_graph().get_each_migration().items()],
                            'origin_node': [[ i, v ] for i, v in sol.node_of_origin.items()]
                        } 
                        for i, sol in enumerate(self)
                    ],
                    'summary': [ [u, v, k]  for (u, v), k in self.summary(draw=False).items()]
                }

    def write_json(self, name='', filename=None):
        if name == '' and filename is not None:
            name = filename.split('.')[0]
        json_dict = self.json(name)
        import json
        if filename is None:
            return json.dumps(json_dict, indent=4)
        else:
            with open(filename, 'w+') as out:
                json.dump(json_dict, out)

    def open_in_viz(self):
        from mach2viz import viz
        viz.Viz(solution=self.json('test')).run()

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