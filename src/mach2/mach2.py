# import os
# import argparse
# import time
# from collections import defaultdict
import gurobipy as gp
from gurobipy import GRB
from collections import defaultdict
from .tree import Refinement
from .solutionset import SolutionSet


class MACH2:
    def __init__(self, clone_tree, primary_location=None, criteria_ordering='MUC', suboptimal_mode=None):

        self.tree = clone_tree
        self.S = clone_tree.locations
        self.V = clone_tree.nodes
        self.E = clone_tree.edges
        self.criteria_ordering = criteria_ordering
        
        self.m = gp.Model('MACH2')

        self._add_vars(criteria_ordering)
        self._add_optimization_function(criteria_ordering)
        self._add_constraints(criteria_ordering)
        self._add_speedup_constraints(criteria_ordering, primary_location)

        if primary_location is not None:
            self._add_primary_location_constraint(primary_location)

    def _add_vars(self, criteria_ordering):
        self.g = self.m.addVars(self.S, self.S, [(v, 'node') for v in self.V] + self.E, vtype=GRB.BINARY, lb=0, ub=1)
        self.r = self.m.addVars(self.S, vtype=GRB.BINARY, lb=0, ub=1)
        self.a = self.m.addVars(self.S, self.S, self.V, vtype=GRB.CONTINUOUS, lb=0, ub=1)
        self.b = self.m.addVars(self.S, self.S, self.V, vtype=GRB.CONTINUOUS, lb=0, ub=1)
        if criteria_ordering[0] == 'C':
            self.z = self.m.addVars(self.S, self.S, vtype=GRB.CONTINUOUS, lb=0, ub=1)

    def _add_optimization_function(self, criteria_ordering):
        sum1 = 0
        for s in self.S:
            for t in self.S:
                if s != t:
                    sum1 += self.g.sum(s, t, '*', '*')
        sum2 = 0
        max_unobserved_clones = 0
        for v in self.V:
            for s in self.S:
                if s not in self.tree.get_labels(v):
                    max_unobserved_clones += 1
                    sum2 += self.g.sum('*', s, '*', v) + self.g.sum('*', s, v, 'node')
                    if v == self.tree.root:
                        sum2 += self.r[s]
        if criteria_ordering in ['MU', 'MUC', 'CMU']:
            self.m.setObjective((max_unobserved_clones + 1) * sum1 + sum2, GRB.MINIMIZE)
        elif criteria_ordering in ['UM', 'UMC', 'CUM']:
            self.m.setObjective((len(self.E) + len(self.V) * (len(self.S) - 1)) * sum2 + sum1, GRB.MINIMIZE)
        elif criteria_ordering in ['M', 'MC', 'CM', 'MCU']:
            self.m.setObjective(sum1, GRB.MINIMIZE)
        elif criteria_ordering in ['U', 'UC', 'CU', 'UCM']:
            self.m.setObjective(sum2, GRB.MINIMIZE)
        else:
            self.m.setObjective(1, GRB.MINIMIZE)

    def _add_constraints(self, criteria_ordering):
        # there is exactly one primary location
        self.m.addConstr( self.r.sum() == 1)

        # for each edge $(u,v)\in E(T)$, there is exactly one edge $(s,t)\in E(R)$ labeled by $e((s,t))=(u,v)$
        for u, v in self.E:
            self.m.addConstr( self.g.sum('*', '*', u, v) == 1 )

        # for each node $u\in V(T)$, $\hat{\ell}(u)\subseteq V(R_u)$.
        for u in self.V:
            for s in self.tree.get_labels(u):
                if u == self.tree.root:
                    self.m.addConstr( self.r[s] + self.g.sum('*', s, u, 'node') == 1 )
                else:
                    self.m.addConstr( self.g.sum('*', s, '*', u) + self.g.sum('*', s, u, 'node') == 1 )

        # each node $t'$ in $R_u$ has at most one incoming edge $(s',t')$ labeled either by $e((s',t'))=u$ or $e((s',t'))=(\pi_T(u),u)$
        for u in self.V:
            for s in self.S:
                if u == self.tree.root:
                    self.m.addConstr( self.r[s] + self.g.sum('*', s, u, 'node') <= 1 )
                else:
                    self.m.addConstr( self.g.sum('*', s, '*', u) + self.g.sum('*', s, u, 'node') <= 1 )

        # set up a and b
        for v in self.V:
            for s in self.S:
                self.m.addConstr(self.a[s, s, v] + self.b[s, s, v] == 0)
                for t in self.S:
                    self.m.addConstr( self.a[s, t, v] >= self.g[s, t, v, 'node'] )
                    self.m.addConstr( self.b[s, t, v] == self.b[t, s, v] )
                    for sp in self.S:
                        self.m.addConstr( self.a[s, t, v] >= self.a[s, sp, v] + self.g[sp, t, v, 'node'] - 1 )
                        if s < t:
                            self.m.addConstr( self.b[s, t, v] >= self.g[sp, s, v, 'node'] + self.g[sp, t, v, 'node'] - 1 )
                            self.m.addConstr( self.b[s, t, v] >= self.b[s, sp, v] + self.g[sp, t, v, 'node'] - 1 )
                    if s < t:
                        s_exists = self.r[s] + self.g.sum('*', s, v, 'node') if v == self.tree.root else self.g.sum('*', s, '*', v) + self.g.sum('*', s, v, 'node')
                        t_exists = self.r[t] + self.g.sum('*', t, v, 'node') if v == self.tree.root else self.g.sum('*', t, '*', v) + self.g.sum('*', t, v, 'node')
                        self.m.addConstr( self.a[s, t, v] + self.a[t, s, v] + self.b[s, t, v] <= s_exists)
                        self.m.addConstr( self.a[s, t, v] + self.a[t, s, v] + self.b[s, t, v] <= t_exists)
                        self.m.addConstr( self.a[s, t, v] + self.a[t, s, v] + self.b[s, t, v] >= s_exists + t_exists - 1)

        # for each node $s'$ in $R_u$, there is a path from $t=r(R_u)$ to $s'$ in $R_u$
        for v in self.V:
            if v == self.tree.root:
                for s in self.S:
                    for t in self.S:
                        self.m.addConstr( self.a[s, t, v] >= self.r[s] + self.g.sum('*', t, v, 'node') - 1 )
            else:
                for s in self.S:
                    for t in self.S:
                        self.m.addConstr( self.a[s, t, v] >= self.g.sum('*', s, '*', v) + self.g.sum('*', t, v, 'node') - 1 )

        for u in self.V:
            if u == self.tree.root:
                for s in self.S:
                    for v in self.tree.get_children(u):
                        self.m.addConstr( self.r[s] + self.g.sum('*', s, u, 'node') >= self.g.sum(s, '*', u, v) )
            else:
                for s in self.S:
                    for v in self.tree.get_children(u):
                        self.m.addConstr( self.g.sum('*', s, '*', u) + self.g.sum('*', s, u, 'node') >= self.g.sum(s, '*', u, v) )

        for v in self.V:
            if v == self.tree.root:
                for s in self.S:
                    if s not in self.tree.get_labels(v):
                        self.m.addConstr( self.g.sum(s, '*', v, '*') >= self.r[s] + self.g.sum('*', s, v, 'node') )
            else:
                for s in self.S:
                    if s not in self.tree.get_labels(v):
                        self.m.addConstr( self.g.sum(s, '*', v, '*') >= self.g.sum('*', s, '*', v) + self.g.sum('*', s, v, 'node') )

        if criteria_ordering[0] == 'C':
            for s in self.S:
                self.m.addConstr( self.z[s,s] == 0)
                self.m.addConstr( self.z.sum('*', s) <= 1 - self.r[s] )
                for t in self.S:
                    if t != s:
                        for v in self.V:
                            self.m.addConstr( self.z[s, t] >= self.g[s, t, v, 'node'] )
                        for u, v in self.E:
                            self.m.addConstr( self.z[s, t] >= self.g[s, t, u, v] )

    def _add_speedup_constraints(self, criteria_ordering, primary_location=None):
        for v in self.V:
            for s in self.S:
                self.m.addConstr( self.g[s, s, v, 'node'] == 0 )
                for t in self.S:
                    self.m.addConstr( self.g[s, t, v, 'node'] + self.g[t, s, v, 'node'] <= 1 )

        if criteria_ordering[0] == 'M':
            for v in self.V:
                sum1 = 0
                for s in self.S:
                    if s not in self.tree.get_labels(v):
                        sum1 += self.g.sum('*', s, '*', v) + self.g.sum('*', s, v, 'node')
                        if v == self.tree.root  and primary_location != s:
                            sum1 += self.r[s]
                if len(self.tree.get_children(v)) < 3:
                    self.m.addConstr( sum1 <= 1 )
                else:
                    self.m.addConstr( sum1 <= len(self.tree.get_children(v)) - 1 )

        if criteria_ordering[0] == 'U':
            for v in self.V:
                if len(self.tree.get_labels(v)) > 0:
                    sum1 = 0
                    for s in self.S:
                        if s not in self.tree.get_labels(v):
                            sum1 += self.g.sum('*', s, '*', v) + self.g.sum('*', s, v, 'node')
                            if v == self.tree.root and primary_location != s:
                                sum1 += self.r[s]
                    self.m.addConstr( sum1 == 0 )
                else:
                    sum1 = 0
                    for s in self.S:
                        sum1 += self.g.sum('*', s, '*', v) + self.g.sum('*', s, v, 'node')
                        if v == self.tree.root:
                            sum1 += self.r[s]
                    self.m.addConstr( sum1 == 1 )


    def _add_primary_location_constraint(self, primary_location):
        self.m.addConstr( self.r[primary_location] == 1 )

    def _count_retrieved_solutions(self):
        self.m.setParam(GRB.Param.SolutionNumber, 0)
        best_obj_val = self.m.PoolObjVal
        for e in range(self.m.SolCount):
            self.m.setParam(GRB.Param.SolutionNumber, e)
            if self.m.PoolObjVal - best_obj_val > 0.5:
                return e
        return self.m.SolCount

    def solve(self, logfile='', starting_nsols = 37, max_solutions=37888, threads=0):
        while True:
            self.m.setParam(GRB.Param.MIPGap, 0)
            self.m.setParam(GRB.Param.PoolSearchMode, 2)
            self.m.setParam(GRB.Param.Threads, threads)
            self.m.setParam(GRB.Param.LogToConsole, 0)
            self.m.setParam(GRB.Param.LogFile, logfile)
            self.m.setParam(GRB.Param.PoolSolutions, starting_nsols)
            self.m.optimize()
            n = self._count_retrieved_solutions()
            if starting_nsols > n:
                self.successfully_run = True
                break
            elif starting_nsols * 2 > max_solutions:
                self.successfully_run = False
                break
            else:
                starting_nsols *= 2
                self.m.reset(1)

        Rs = []
        for e in range(n):
            R = defaultdict(list)
            self.m.setParam(GRB.Param.SolutionNumber, e)
            for s in self.S:
                for t in self.S:
                    for v in self.V:
                        if self.g[s, t, v, 'node'].Xn > 0.5:
                            R[v].append((s,t))
                    for u, v in self.E:
                        if self.g[s, t, u, v].Xn > 0.5:
                            R[(u,v)].append((s,t))
            Rs.append(R)
        
        return SolutionSet([Refinement.from_refinement_graph(self.tree, R) for R in Rs]).filter(criteria_ordering=self.criteria_ordering)
