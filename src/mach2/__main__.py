import argparse
import os
import time
from  .tree import MultiLabeledTree, Refinement
from .mach2 import MACH2
from . import utils

def process_args():
    parser = argparse.ArgumentParser(description='MACH2')

    parser.add_argument('clonal_tree', type=str, help='Input clonal tree')
    parser.add_argument('observed_labeling', type=str, help='Input observed labeling')

    parser.add_argument('-c', '--criteria', type=str, help='Criteria ordering', default='UMC')
    parser.add_argument('-s', '--seeding_locations', action='store_true', default=False, 
                help='Prioritize solutions with the least number of seeding locations (default=False)')
    parser.add_argument('-p', '--primary', type=str, help='Primary anatomical location')
    parser.add_argument('--colormap', metavar='COLORMAP', type=str, help='Color map file', action='store')
    parser.add_argument('--log', action='store_true', default=False, help='Outputs Gurobi logging (default=False)')
    parser.add_argument('-o', '--output', action='store', default=None, help='Output folder (default=current folder)')
    parser.add_argument('--max_solutions', type=int, help='Maximum number of solutions retained (default=37888)', default=37888)
    parser.add_argument('-t', '--threads', type=int, help='Number of threads', default=0)
    # parser.add_argument('-s', '--suboptimal', action='store_true', default=False, help='Returns suboptimal solutions without duplicates, \
        # may be slow (default=False)')
    parser.add_argument('--viz', '--open_in_viz', action='store_true', default=False, help='Open the locations on MACH2-viz \
        (default=False)')

    return parser.parse_args()

def main():
    args = process_args()
    tree = MultiLabeledTree.from_files(args.clonal_tree, args.observed_labeling)

    if args.output is None:
        output_str = '.'
    else:
        output_str = args.output
        if not os.path.exists(output_str):
            os.makedirs(output_str)

    if args.primary is None:
        primary_str = 'ALL'
    else:
        primary_str = args.primary

    if args.log == True:
        logfile = f'{output_str}/{primary_str}-log.txt'
    else:
        logfile = ''

    start_t = time.time()
    solver = MACH2(tree, primary_location=args.primary, criteria_ordering=args.criteria)
    solutions = solver.solve(logfile=logfile,  max_solutions=args.max_solutions, threads=args.threads)
    if args.seeding_locations:
        solutions = solutions.seeding_location_filter()
    total_t = time.time() - start_t

    if args.viz:
        solutions.open_in_viz()
    else:
        if args.colormap:
            colormap = utils.process_colormap_file(args.colormap)
        else:
            colormap = utils.get_colormap(phylogeny.locations)

        padding = len(str(len(solutions)))

        for e, soln in enumerate(solutions):
            primary_str = soln.get_label(soln.root)
            soln.write_refinement(f'{output_str}/{primary_str}-T-{str(e).zfill(padding)}')
            soln.migration_graph().write_graph(f'{output_str}/{primary_str}-G-{str(e).zfill(padding)}.graph')
            soln.draw(f'{output_str}/{primary_str}-T-{str(e).zfill(padding)}.dot')
            soln.migration_graph().draw(f'{output_str}/{primary_str}-G-{str(e).zfill(padding)}.dot')
            print(f'{primary_str}-\t{e}\t{len(soln.unobserved_clones)}\t{len(soln.migrations)}\t{len(soln.comigrations)}\t{total_t}')
        