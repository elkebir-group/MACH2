import networkx as nx
import pandas as pd
from collections import defaultdict
import glob
from .tree import Refinement
import graphviz as gv


class SolutionSet:
    """
    A class representing a set of solutions (Refinement objects) derived from a common unrefined tree.
    Provides methods for filtering, summarizing, and visualizing the solutions.
    """
    def __init__(self, sol_list=[], check=True, ranked=True):
        """
        Initializes the SolutionSet with a list of Refinement objects.
        
        :param sol_list: A list of Refinement objects.
        :param check: If True, ensures all solutions share the same unrefined tree.
        :param ranked: If True, maintains the order of solutions as provided.
        """
        self.solution_set = set(sol_list)
        self.unrefined_tree = [i for i in sol_list][0].unrefined_tree if len(sol_list) > 0 else None
        if ranked:
            self.sol_list = sol_list
        if check:
            for sol in self.solution_set:
                if self.unrefined_tree != sol.unrefined_tree:
                    raise ValueError('Solutions with different original unrefined trees')
                
    def from_files(output_prefix, utreefile, ulabelfile, colormap_file=None):
        """
        Creates a SolutionSet from files containing refined trees and associated data.
        
        :param output_prefix: Prefix for the output files.
        :param utreefile: Path to the unrefined tree file.
        :param ulabelfile: Path to the unrefined label file.
        :param colormap_file: Optional path to the colormap file.
        :return: A SolutionSet object.
        """
        sol_list = []
        for treefile in glob.glob(output_prefix + '*.refined.tree'):
            sol_list.append(Refinement.from_files(treefile, 
                                                  treefile[:-13]+'.location.labeling',
                                                  utreefile, ulabelfile,
                                                  treefile[:-13]+'.nodes',
                                                  colormap_file))
        return SolutionSet(sol_list, check=False)

    def __len__(self):
        """Returns the number of solutions in the set."""
        return len(self.solution_set)
    
    def __iter__(self):
        """Returns an iterator over the solutions in the set."""
        return iter(self.solution_set)

    def __getitem__(self, index):
        """Returns the solution at the specified index."""
        if not hasattr(self, 'sol_list'):
            self.sol_list = list(self.solution_set)
        return self.sol_list[index]
    
    def __add__(s1, s2):
        """
        Combines two SolutionSets into one.
        
        :param s1: The first SolutionSet.
        :param s2: The second SolutionSet.
        :return: A new SolutionSet containing all solutions from both sets.
        """
        if len(s1.solution_set) == 0:
            return s2
        elif len(s2.solution_set) == 0:
            return s1
        return SolutionSet(s1.solution_set.union(s2.solution_set), check=(s1.unrefined_tree!=s2.unrefined_tree))

    def __sub__(s1, s2):
        """
        Subtracts one SolutionSet from another.
        
        :param s1: The first SolutionSet.
        :param s2: The second SolutionSet.
        :return: A new SolutionSet containing solutions in s1 but not in s2.
        """
        if len(s1.solution_set) == 0:
            return SolutionSet()
        elif len(s2.solution_set) == 0:
            return s1
        return SolutionSet(s1.solution_set - s2.solution_set, check=(s1.unrefined_tree!=s2.unrefined_tree), ranked=False)

    def intersection(s1, s2):
        """
        Returns the intersection of two SolutionSets.
        
        :param s1: The first SolutionSet.
        :param s2: The second SolutionSet.
        :return: A new SolutionSet containing solutions common to both sets.
        """
        if len(s1.solution_set) == 0 or len(s2.solution_set) == 0:
            return SolutionSet()
        return SolutionSet(s1.solution_set.intersection(s2.solution_set), check=(s1.unrefined_tree!=s2.unrefined_tree))
    
    def filter(self, criteria_ordering):
        """
        Filters the solutions based on the specified criteria ordering.
        
        :param criteria_ordering: A string specifying the order of optimization criteria (e.g., 'UMCS').
        :return: A new SolutionSet containing only the solutions that match the criteria.
        """
        if criteria_ordering is not None:
            sort_key_map = { 
                'U': lambda ref: len(ref.unobserved_clones), 
                'M': lambda ref: len(ref.migrations), 
                'C': lambda ref: len(ref.comigrations),
                'S': lambda ref: len(ref.seeding_locations),
            }
            sorter = lambda ref, co: tuple(sort_key_map[letter](ref) for letter in co)
            min_val = min([sorter(sol, criteria_ordering) for sol in self.solution_set])
            return SolutionSet([sol for sol in self.solution_set if sorter(sol, criteria_ordering) == min_val], check=False)
        else:
            return self
        

    def co_occurence_table(self):
        """
        Generates a co-occurrence table of migrations across all solutions.
        
        :return: A pandas DataFrame representing the co-occurrence table.
        """
        table = defaultdict(lambda: defaultdict(int))
        for solution in self:
            for u1, v1 in solution.migration_graph().get_migrations():
                for u2, v2 in solution.migration_graph().get_migrations():
                    table[f"{u1} -> {v1}"][f"{u2} -> {v2}"] += 1
        return pd.DataFrame(table).fillna(0)
    
    def summary(self, draw=True):
        """
        Generates a summary of migrations across all solutions.
        
        :param draw: If True, returns a Graphviz visualization of the summary.
        :return: A Graphviz Digraph object if draw=True, otherwise a dictionary of migration counts.
        """
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
        """
        Writes all solutions in the set to files.
        
        :param output_prefix: Prefix for the output files.
        """
        for i, sol in enumerate(self):
            sol.write_refinement(f'{output_prefix}sol{i}')

    def json(self, name=''):
        """
        Converts the SolutionSet to a JSON-compatible dictionary. Compatible with MACH2-viz.
        
        :param name: Optional name for the solution set.
        :return: A dictionary representing the solution set in JSON format.
        """
        return {
                    'name' : name,
                    'original' : {
                        'tree': [[u,v] for u,v in self.unrefined_tree.edges], 
                        'labeling': [ [v, list(self.unrefined_tree.get_labels(v))] for v in self.unrefined_tree.nodes]
                    },
                    'solutions' : [
                        {
                            'name' : f'\\#{i}',
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
        """
        Writes the SolutionSet to a JSON file or returns the JSON string.
        
        :param name: Optional name for the solution set.
        :param filename: Optional path to the output JSON file.
        :return: If filename is None, returns the JSON string. Otherwise, writes to the file.
        """
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
        """
        Opens the SolutionSet in the MACH2-viz visualization tool.
        
        :raises ImportError: If MACH2-viz is not installed.
        """
        try:
            from mach2viz import viz
            viz.Viz(solution=self.json('Direct_from_MACH2')).run()
        except:
            raise ImportError('Install MACH2-viz to use this feature.')