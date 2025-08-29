try:
    from .mach2 import MACH2
except ModuleNotFoundError:
    print('Gurobi not found. MACH2 requires Gurobi to run.\nYou can still analyze clonal trees and migration graphs.')
from  .tree import MultiLabeledTree, Refinement
from .solutionset import SolutionSet