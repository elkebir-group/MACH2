import pyomo.environ as pyo
from pyomo.opt import SolverFactory
import networkx as nx
import os
import argparse
import time
from collections import defaultdict
from .clonetree import CloneTree, RefinedCloneTree
from .migrationgraph import MigrationGraph
from .solutionset import SolutionSet
from . import utils


class MACH:

    def __init__(self, clone_tree, primary_site = None, suboptimal_mode=False, seeding_site=False, specify_migration_comigration=None, possible_migration_list=None):
        if specify_migration_comigration is not None:
            suboptimal_mode = True
        self.clone_tree = clone_tree
        self.E = self.clone_tree.edges
        self.V = [(i, 'node') for i in self.clone_tree.nodes]
        self.L = [(i, 'node') for i in self.clone_tree.leaves]
        self.Sigma = self.clone_tree.sites
        self.suboptimal_mode = suboptimal_mode
        self.seeding_site = seeding_site
        self.specific_mig_comig = specify_migration_comigration

        self.m = pyo.ConcreteModel()

        self._add_vars()
        if suboptimal_mode:
            self._add_extra_vars_suboptimal()
        if seeding_site:
            self._add_extra_vars_seeding_sites()

        if not seeding_site:
            self._add_optimization_function()
        else:
            self._add_optimization_function_seeding_sites()

        self._add_leaf_constraints()
        self._add_root_constraints()
        self._add_constraints_for_p()
        self.add_polytomy_resolution_compatibility_constraints()
        self._add_original_edges_compatibility_constraints()
        clone_tree.infer_paths()
        if suboptimal_mode:
            self.add_constraints_for_z_b_R_suboptimal()
            if seeding_site:
                self.add_constraints_for_q_suboptimal()
        else:
            self.add_constraints_for_z()
            if seeding_site:
                self.add_constraints_for_q()
        self.add_constraints_binary()
        if primary_site is not None:
            self.add_primary_site_constraints(primary_site)
        if specify_migration_comigration is not None:
            self.mig, self.comig = specify_migration_comigration
            self.add_contraints_specifying_mig_comig()
        if possible_migration_list is not None:
            self._add_constraints_2_pick_migrations_4m_given_list(possible_migration_list)

    def _add_vars(self):
        self.m.l = pyo.Var(self.V, self.Sigma, domain=pyo.Binary)
        self.m.g = pyo.Var(self.Sigma, self.Sigma,
                           self.E + self.V, domain=pyo.Binary)
        self.m.z = pyo.Var(self.Sigma, self.Sigma,
                           domain=pyo.NonNegativeIntegers)
        self.m.r = pyo.Var(self.Sigma, domain=pyo.Binary)
        self.m.p = pyo.Var(self.Sigma, self.Sigma, self.E, domain=pyo.Binary)

    def _add_extra_vars_suboptimal(self):
        self.m.R = pyo.Var(self.Sigma, self.Sigma, self.L, domain=pyo.Binary)
        self.m.b = pyo.Var(self.Sigma, self.Sigma, self.L, domain=pyo.Binary)

    def _add_extra_vars_seeding_sites(self):
        self.m.q = pyo.Var(self.Sigma, domain=pyo.Binary)

    def _add_optimization_function(self):
        self.m.obj = pyo.Objective(expr=sum(self.m.g[s, t, k] for s in self.Sigma for t in self.Sigma for k in self.E + self.V if s != t) +
                                   sum(self.m.z[s, t] for s in self.Sigma for t in self.Sigma if s != t))

    def _add_optimization_function_seeding_sites(self):
        self.m.obj = pyo.Objective(expr=sum(self.m.g[s, t] for s in self.Sigma for t in self.Sigma if s != t) +
                                   sum(self.m.z[s, t] for s in self.Sigma for t in self.Sigma if s != t) +
                                   sum(self.m.q[s] for s in self.Sigma))

    def _add_leaf_constraints(self):
        self.m.leaf_constraints = pyo.ConstraintList()
        for i in self.L:
            self.m.leaf_constraints.add(
                self.m.l[i, self.clone_tree.get_label(i[0])] == 1)
            self.m.leaf_constraints.add(sum(
                self.m.l[i, s] for s in self.Sigma if s != self.clone_tree.get_label(i[0])) == 0)
            
    def _add_root_constraints(self):
        self.m.root_constraints = pyo.ConstraintList()
        self.m.root_constraints.add(sum(self.m.r[s] for s in self.Sigma) == 1)

    def _add_constraints_for_p(self):
        self.m.p_constraints = pyo.ConstraintList()
        for s in self.Sigma:
            for t in self.Sigma:
                for i, j in self.E:
                    sum1 = sum(self.m.g[t, t_prime, (i, j)]
                               for t_prime in self.Sigma)
                    sum2 = sum(self.m.p[t, t_prime, (i, j)]
                               for t_prime in self.Sigma)

                    self.m.p_constraints.add(
                        self.m.p[s, t, (i, j)] >= self.m.g[s, t, (i, 'node')] + sum1 - 1)
                    self.m.p_constraints.add(
                        self.m.p[s, t, (i, j)] >= self.m.g[s, t, (i, 'node')] + sum2 - 1)
                    self.m.p_constraints.add(
                        self.m.p[s, t, (i, j)] <= self.m.g[s, t, (i, 'node')])
                    self.m.p_constraints.add(
                        self.m.p[s, t, (i, j)] <= sum1 + sum2)

    def add_polytomy_resolution_compatibility_constraints(self):

        self.m.polytomy_constraints_set_l = pyo.ConstraintList()

        for i in self.V:
            for s in self.Sigma:
                for t in self.Sigma:
                    if s != t:
                        self.m.polytomy_constraints_set_l.add(
                            self.m.g[s, t, i] <= self.m.l[i, s])
                        self.m.polytomy_constraints_set_l.add(
                            self.m.g[s, t, i] <= self.m.l[i, t])
                    else:
                        self.m.polytomy_constraints_set_l.add(
                            self.m.g[s, t, i] == 0)

        for i in self.V:
            for s in self.Sigma:
                sum1 = sum(self.m.g[t, s, i] for t in self.Sigma)
                if i[0] != self.clone_tree.root:
                    pii = self.clone_tree.get_parent_arc(i[0])
                    sum1 += sum(self.m.g[t, s, pii] for t in self.Sigma)
                else:
                    sum1 += self.m.r[s]
                self.m.polytomy_constraints_set_l.add(sum1 == self.m.l[i, s])

        for i in self.V:
            delta_i = self.clone_tree.get_children_arcs(i[0])
            for s in self.Sigma:
                for t in self.Sigma:
                    if s != t:
                        self.m.polytomy_constraints_set_l.add(self.m.g[s, t, i] <= sum(self.m.p[s, t, ij] for ij in delta_i))

        for i in self.V:
            if i[0] != self.clone_tree.root:
                pii = self.clone_tree.get_parent_arc(i[0])
            delta_i = self.clone_tree.get_children_arcs(i[0])
            for ij in delta_i:
                for s in self.Sigma:
                    sum1 = 0
                    if i[0] != self.clone_tree.root:
                        sum1 = sum(self.m.g[t,s,pii] for t in self.Sigma)
                    else:
                        sum1 = self.m.r[s]
                    sum2 = sum(self.m.p[s,t,ij]+self.m.g[s,t,ij] for t in self.Sigma)
                    self.m.polytomy_constraints_set_l.add(sum1 <= sum2)

        # self.m.polytomy_constraints_valid_structure = pyo.ConstraintList()
        # for i in self.V:
        #     sum1 = sum(self.m.g[s, t, i]
        #                for s in self.Sigma for t in self.Sigma)
        #     sum2 = sum(self.m.l[i, s] for s in self.Sigma)
        #     self.m.polytomy_constraints_valid_structure.add(sum1 == sum2 - 1)

        # for s in self.Sigma:
        #     for ij in self.E:
        #         self.m.polytomy_constraints_valid_structure.add(
        #             sum(self.m.g[s, t, ij] + self.m.p[s, t, ij] for t in self.Sigma) <= 1)

        # for i in self.V:
        #     for s in self.Sigma:
        #         if i[0] != self.clone_tree.root:
        #             pii = self.clone_tree.get_parent_arc(i[0])
        #         delta_i = self.clone_tree.get_children_arcs(i[0])
        #         sum1 = 0
        #         for t in self.Sigma:
        #             if i[0] != self.clone_tree.root:
        #                 sum1 += self.m.g[t, s, pii]
        #             for ij in delta_i:
        #                 sum1 += self.m.g[s, t, ij]
        #         self.m.polytomy_constraints_valid_structure.add(
        #             sum1 >= self.m.l[i, s])

        for i in self.V:
            for s in self.Sigma:
                sum1 = 0
                sum2 = 0
                for t in self.Sigma:
                    if s != t:
                        sum1 += self.m.g[s, t, i] + self.m.g[t, s, i]
                        sum2 += self.m.l[i, t]
                self.m.polytomy_constraints_set_l.add(
                    self.clone_tree.n_sites * sum1 >= sum2 + self.clone_tree.n_sites * self.m.l[i, s] - self.clone_tree.n_sites)

    def _add_original_edges_compatibility_constraints(self):
        self.m.orig_edge = pyo.ConstraintList()
        for s in self.Sigma:
            for t in self.Sigma:
                for (i, j) in self.E:
                    self.m.orig_edge.add(
                        self.m.g[s, t, (i, j)] <= self.m.l[(i, 'node'), s])
                    self.m.orig_edge.add(
                        self.m.g[s, t, (i, j)] <= self.m.l[(j, 'node'), t])

        for ij in self.E:
            sum1 = 0
            for s in self.Sigma:
                for t in self.Sigma:
                    sum1 += self.m.g[s, t, ij]
            self.m.orig_edge.add(sum1 == 1)

    def add_constraints_for_z(self):
        self.m.constraints_z = pyo.ConstraintList()
        for s in self.Sigma:
            for t in self.Sigma:
                if s != t:
                    for k in self.L:
                        sum1 = 0
                        path_k = self.clone_tree.paths[k[0]]
                        for uv in path_k:
                            sum1 += self.m.g[s, t, uv]
                            sum1 += self.m.p[s, t, uv]
                        self.m.constraints_z.add(self.m.z[s, t] >= sum1)
                else:
                    self.m.constraints_z.add(self.m.z[s, t] == 0)

    def add_constraints_for_z_b_R_suboptimal(self):
        self.m.constraints_zbr = pyo.ConstraintList()
        for s in self.Sigma:
            for t in self.Sigma:
                if s != t:
                    sum3 = 0
                    for k in self.L:
                        path_k = self.clone_tree.paths[k[0]]
                        sum1 = sum((self.m.g[s, t, uv] + self.m.p[s, t, uv]) for uv in path_k)
                        self.m.constraints_zbr.add(self.m.z[s, t] >= sum1)
                        self.m.constraints_zbr.add(
                            self.m.z[s, t] <= sum1 + self.clone_tree.max_height * (1 - self.m.b[s, t, k]))
                        self.m.constraints_zbr.add(
                            self.m.z[s, t] - sum1 >= 1 - self.m.R[s, t, k])
                        self.m.constraints_zbr.add(
                            self.m.z[s, t] - sum1 <= self.clone_tree.max_height * (1 - self.m.R[s, t, k]))
                        self.m.constraints_zbr.add(
                            self.m.b[s, t, k] >= self.m.R[s, t, k] - sum3)
                        sum3 += self.m.R[s, t, k]
                    sum2 = sum(self.m.b[s, t, k] for k in self.L)
                    self.m.constraints_zbr.add(sum2 == 1)
                else:
                    self.m.constraints_zbr.add(self.m.z[s, t] == 0)
                    for k in self.L:
                        self.m.constraints_zbr.add(self.m.R[s, t, k] == 0)
                        self.m.constraints_zbr.add(self.m.b[s, t, k] == 0)

    def add_constraints_for_q(self):
        self.m.constraints_s = pyo.ConstraintList()
        for s in self.Sigma:
            self.m.constraints_s.add(self.m.q[s] <= 1)
            for t in self.Sigma:
                if s != t:
                    for i in self.V + self.E:
                        self.m.constraints_s.add(
                            self.m.q[s] >= self.m.g[s, t, i])

    def add_constraints_for_q_suboptimal(self):
        self.m.constraints_s_suboptimal = pyo.ConstraintList()
        for s in self.Sigma:
            self.m.constraints_s_suboptimal.add(self.m.q[s] <= 1)
            sum1 = 0
            for t in self.Sigma:
                if s != t:
                    for i in self.V + self.E:
                        self.m.addConstr(self.m.q[s] >= self.m.g[s, t, i])
                        sum1 += self.m.g[s, t, i]
            self.m.constraints_s_suboptimal.add(self.q[s] <= sum1)

    def add_constraints_binary(self):
        self.m.constraints_binary = pyo.ConstraintList()
        for i in self.V:
            delta_i = self.clone_tree.get_children_arcs(i[0])
            if len(delta_i) <= 2:
                sum1 = 0
                for s in self.Sigma:
                    sum1 += self.m.l[i, s]
                self.m.constraints_binary.add(sum1 == 1)
            else:
                for s in self.Sigma:
                    sum1 = 0
                    for t in self.Sigma:
                        sum1 += self.m.g[s, t, i]
                        for ij in delta_i:
                            sum1 += self.m.g[s, t, ij]
                    self.m.constraints_binary.add(sum1 >= 2 * self.m.l[i, s])

    def add_primary_site_constraints(self, primary_site):
        self.m.constraints_primary = pyo.ConstraintList()
        self.m.constraints_primary.add(self.m.r[primary_site] == 1)

    def add_contraints_specifying_mig_comig(self):
        self.m.constraints_mig_comig = pyo.ConstraintList()
        self.m.constraints_mig_comig.add(sum(
            self.m.g[s, t, k] for s in self.Sigma for t in self.Sigma for k in self.E + self.V if s != t) == self.mig)
        self.m.constraints_mig_comig.add(sum(
            self.m.z[s, t] for s in self.Sigma for t in self.Sigma if s != t) == self.comig)
        
    def _add_constraints_2_pick_migrations_4m_given_list(self, possible_migration_list):
        self.m.constraints_mig_list = pyo.ConstraintList()
        for s in self.Sigma:
            for t in self.Sigma:
                if (s, t) not in possible_migration_list:
                    self.m.constraints_mig_list.add(self.m.z[s, t] == 0)
        
    def _count_optimal_solution(self, opt):
        opt._solver_model.params.SolutionNumber = 0
        best_obj_val = opt._solver_model.PoolObjVal
        for e in range(opt._solver_model.SolCount):
            opt._solver_model.params.SolutionNumber = e
            if opt._solver_model.PoolObjVal - best_obj_val > 0.5:
                return e
        return opt._solver_model.SolCount

    def solve(self, solver, nSolutions, logfile=None, n_threads=1, raw=False):
        if solver == 'gurobi':
            opt = SolverFactory("gurobi", solver_io='python', options={ 'MIPGap': 0, 'PoolSolutions': nSolutions, 'PoolSearchMode': 2, 'threads': n_threads, 'LogToConsole': 0, 'LogFile': logfile})
            opt.solve(self.m, load_solutions=True, tee=True)
            if self.suboptimal_mode and self.specific_mig_comig is None:
                n_actual_solutions = opt._solver_model.SolCount
            else:
                n_actual_solutions = self._count_optimal_solution(opt)

            solutions_raw = []
            for soln in range(n_actual_solutions):
                opt._solver_model.params.SolutionNumber = soln
                ell = defaultdict(list)
                G = defaultdict(list)
                r = []
                z = []
                p = defaultdict(list)
                b = defaultdict(list)
                mu = 0
                gamma = 0
                for i in self.V:
                    for s in self.Sigma:
                        if opt._pyomo_var_to_solver_var_map[self.m.l[i, s]].Xn > 0.5:
                            ell[i[0]].append(s)
                        # print(f'l[{i},{s}] = {opt._pyomo_var_to_solver_var_map[self.m.l[i, s]].Xn}')
                for s_1 in self.Sigma:
                    for s_2 in self.Sigma:
                        for i in self.E:
                            if opt._pyomo_var_to_solver_var_map[self.m.g[s_1, s_2, i]].Xn > 0.5:
                                G[i].append((s_1, s_2))
                                if s_1 != s_2:
                                    mu += 1
                            if opt._pyomo_var_to_solver_var_map[self.m.p[s_1, s_2, i]].Xn > 0.5:
                                p[i].append((s_1, s_2))
                            # print(f'G[{s_1},{s_2},{i}] = {opt._pyomo_var_to_solver_var_map[self.m.g[s_1, s_2, i]].Xn}')
                            # print(f'p[{s_1},{s_2},{i}] = {opt._pyomo_var_to_solver_var_map[self.m.p[s_1, s_2, i]].Xn}')
                        for i in self.V:
                            if opt._pyomo_var_to_solver_var_map[self.m.g[s_1, s_2, i]].Xn > 0.5:
                                G[i[0]].append((s_1, s_2))
                                if s_1 != s_2:
                                    mu += 1
                            # print(f'G[{s_1},{s_2},{i}] = {opt._pyomo_var_to_solver_var_map[self.m.g[s_1, s_2, i]].Xn}')
                for s_1 in self.Sigma:
                    for s_2 in self.Sigma:
                        if opt._pyomo_var_to_solver_var_map[self.m.z[s_1, s_2]].Xn > 0.5:
                            z.append((s_1,s_2, opt._pyomo_var_to_solver_var_map[self.m.z[s_1, s_2]].Xn))
                            gamma += 1#opt._pyomo_var_to_solver_var_map[self.m.z[s_1, s_2]].Xn
                        # for l in self.L:
                            # if opt._pyomo_var_to_solver_var_map[self.m.b[s_1, s_2, l]].Xn > 0.5:
                                # b[(s_1,s_2)].append(l)
                        # print(f'z[{s_1},{s_2}] = {opt._pyomo_var_to_solver_var_map[self.m.z[s_1, s_2]].Xn}')
                for s in self.Sigma:
                    if opt._pyomo_var_to_solver_var_map[self.m.r[s]].Xn > 0.5:
                        r.append(s)
                    # print(f'r[{s}] = {opt._pyomo_var_to_solver_var_map[self.m.r[s]].Xn}')
                solution_raw = {'vertex_multilabeling': ell,
                                'aug_mig_graph': G, 'n_mig': int(mu), 'n_comig': int(gamma),#}
                                  'other': {'r': r, 'z': z, 'p':p, 'b':b}}
                if self.seeding_site:
                    sigma = 0
                    for s in self.Sigma:
                        if opt._pyomo_var_to_solver_var_map[self.m.q[s]].Xn > 0.5:
                            sigma += 1
                    solution_raw['n_seedingsites'] = sigma
                solutions_raw.append(solution_raw)
        if raw:
            return solutions_raw
        else:
            return SolutionSet([(RefinedCloneTree(self.clone_tree, raw), MigrationGraph(raw)) for raw in solutions_raw])


def process_args():
    parser = argparse.ArgumentParser(description='MACH2')

    parser.add_argument('clone_tree', type=str, help='Input clone tree')
    parser.add_argument('leaf_labeling', type=str, help='Input leaf labeling')

    parser.add_argument('-p', '--primary', type=str, help='Primary anatomical site')
    parser.add_argument('-c', '--colormap', metavar='COLORMAP', type=str, help='Color map file', action='store')
    parser.add_argument('--log', action='store_true', default=False, help='Outputs Gurobi logging')
    parser.add_argument('-o', '--output', action='store', default=None, help='Output folder')
    parser.add_argument('-N', '--nsolutions', type=int, help='Maximum number of solutions retained', default=10)
    parser.add_argument('-C', '--count_solutions', action='store_true', default=False, help='Only prints the number of solutions\
        (default=False)')
    parser.add_argument('-t', '--threads', type=int, help='Number of threads')
    parser.add_argument('-s', '--suboptimal', action='store_true', default=False, help='Returns suboptimal solutions without duplicates, \
        may be slow (default=False)')
    parser.add_argument('-S', '--seeding_sites', action='store_true', default=False, help='Minimizes the number of seeding sites \
        too (default=False)')

    return parser.parse_args()

def main():
    args = process_args()
    clone_tree = CloneTree.from_file(args.clone_tree, args.leaf_labeling)

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
    solver = MACH(clone_tree, primary_site=args.primary, suboptimal_mode=args.suboptimal, seeding_site=args.seeding_sites)
    solutions = solver.solve('gurobi', args.nsolutions, logfile=logfile, n_threads=args.threads, raw=False)
    total_t = time.time() - start_t

    if args.count_solutions:
        if args.seeding_sites:
            print(f'{primary_str}-\t{int(solutions[0].n_migrations)}\t{int(solutions[0].n_comigrations)}\t{int(solutions[0].n_seeding_sites)}\t\
                {len(solutions)}')
        else:
            print(f'{primary_str}-\t{int(solutions[0].n_migrations)}\t{int(solutions[0].n_comigrations)}\t{len(solutions)}')
    else:
        if args.colormap:
            colormap = utils.process_colormap(colormap_filename=args.colormap)
        else:
            colormap = utils.get_colormap(clone_tree.sites)
        
        padding = len(str(len(solutions)))

        for e, soln in enumerate(solutions):
            primary_str = soln.clone_tree.primary_site
            soln.clone_tree.write_dot(f'{output_str}/{primary_str}-T-{str(e).zfill(padding)}.dot', colormap=colormap)
            soln.clone_tree.write_tree(f'{output_str}/{primary_str}-T-{str(e).zfill(padding)}.tree')
            soln.clone_tree.write_labeling(f'{output_str}/{primary_str}-T-{str(e).zfill(padding)}.labeling')
            soln.migration_graph.write_dot(f'{output_str}/{primary_str}-G-{str(e).zfill(padding)}.dot', colormap=colormap)
            soln.migration_graph.write_graph(f'{output_str}/{primary_str}-G-{str(e).zfill(padding)}.graph')
            print(f'{primary_str}-\t{e}\t{soln.n_migrations}\t{soln.n_comigrations}\t{total_t}')

if __name__ == '__main__':
    main()