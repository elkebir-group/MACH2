import gurobipy as gp
from gurobipy import GRB
from collections import defaultdict
from .tree import Refinement
from .solutionset import SolutionSet


class MACH2:
    """
    A class representing the MACH2 algorithm for inferring migration histories of metastatic cancer.
    """
    def __init__(self, clone_tree, primary_location=None, criteria_ordering='UMCS'):
        """
        Initializes the MACH2 solver with a clone tree, optional primary location, and criteria ordering.
        
        :param clone_tree: The input clone tree to be refined.
        :param primary_location: Optional primary location constraint.
        :param criteria_ordering: String specifying the order of optimization criteria (e.g., 'UMCS').
        """
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
        """Private function for internal use by the MACH2 algorithm."""
        self.g = self.m.addVars(self.S, self.S, [(v, 'node') for v in self.V] + self.E, vtype=GRB.BINARY, lb=0, ub=1)
        self.r = self.m.addVars(self.S, vtype=GRB.BINARY, lb=0, ub=1)
        self.a = self.m.addVars(self.S, self.S, self.V, vtype=GRB.CONTINUOUS, lb=0, ub=1)
        self.b = self.m.addVars(self.S, self.S, self.V, vtype=GRB.CONTINUOUS, lb=0, ub=1)
        if 'C' in criteria_ordering:
            self.c = self.m.addVars(self.S, self.S, self.E, vtype=GRB.CONTINUOUS, lb=0, ub=1)
            self.d = self.m.addVars(self.S, self.S, [(s, u) for u in self.V for s in self.tree.get_labels(u)], vtype=GRB.CONTINUOUS, lb=0, ub=1)
            self.z = self.m.addVars(self.S, self.S, vtype=GRB.CONTINUOUS, lb=0)
        if 'S' in criteria_ordering:
            self.s = self.m.addVars(self.S, vtype=GRB.BINARY, lb=0, ub=1)

    def _add_optimization_function(self, criteria_ordering):
        """Private function for internal use by the MACH2 algorithm."""
        pscores = {}
        if 'U' in criteria_ordering:
            pscores['U'] = 0
            max_unobserved_clones = 0
            for v in self.V:
                for s in self.S:
                    if s not in self.tree.get_labels(v):
                        max_unobserved_clones += 1
                        pscores['U'] += self.g.sum('*', s, '*', v) + self.g.sum('*', s, v, 'node')
                        if v == self.tree.root:
                            pscores['U'] += self.r[s]
        if 'M' in criteria_ordering:
            pscores['M'] = 0
            for s in self.S:
                for t in self.S:
                    if s != t:
                        pscores['M'] += self.g.sum(s, t, '*', '*')
        if 'C' in criteria_ordering:
            pscores['C'] = self.z.sum()
        if 'S' in criteria_ordering:
            pscores['S'] = self.s.sum()
        score = 0
        for c in criteria_ordering:
            if c == 'U':
                score = score * (max_unobserved_clones + 1) + pscores[c]
            elif c == 'M' or c == 'C':
                score = score * (len(self.E) + len(self.V) * (len(self.S) - 1)) + pscores[c]
            elif c == 'S':
                score = score * len(self.S) + pscores[c]
        self.m.setObjective(score, GRB.MINIMIZE)

    def _add_constraints(self, criteria_ordering):
        """Private function for internal use by the MACH2 algorithm."""
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
        
        if 'C' in criteria_ordering:
            for s in self.S:
                for t in self.S:
                    for (u,v) in self.E:
                        self.m.addConstr( self.c[s, t, u, v] >= self.g[s, t, u, v])
                        self.m.addConstr( self.c[s, t, u, v] >= self.g[s, t, u, 'node'] + self.g.sum(t, '*', u, v)  - 1 )
                        self.m.addConstr( self.c[s, t, u, v] <= self.g[s, t, u, v] + self.g[s, t, u, 'node'] )
                        self.m.addConstr( self.c[s, t, u, v] <= self.g[s, t, u, v] + self.g.sum(t, '*', u, v) )
                    for u in self.V:
                        sum1 = 0
                        for (u1, v1) in self.tree.paths[u]:
                            sum1 += self.c[s, t, u1, v1]
                        for sp in self.tree.get_labels(u):
                            if t == sp:
                                self.m.addConstr( self.d[s, t, sp, u] == self.g[s, t, u, 'node'] )
                            else:
                                self.m.addConstr( self.d[s, t, sp, u] <= self.g[s, t, u, 'node'] )
                                self.m.addConstr( self.d[s, t, sp, u] <= self.a[t, sp, u] )
                                self.m.addConstr( self.d[s, t, sp, u] >= self.g[s, t, u, 'node'] + self.a[t, sp, u] - 1 )
                            if s == t:
                                self.m.addConstr( self.z[s, t] == 0 )
                            else:
                                self.m.addConstr( self.z[s, t] >= sum1 + self.d[s, t, sp, u] )
        
        if 'S' in criteria_ordering:
            for u, v in [(v, 'node') for v in self.V] + self.E:
                for s in self.S:
                    for t in self.S:
                        if s != t:
                            self.m.addConstr( self.s[s] >= self.g[s, t, u, v] )


    def _add_speedup_constraints(self, criteria_ordering, primary_location=None):
        """Private function for internal use by the MACH2 algorithm."""
        for v in self.V:
            for s in self.S:
                self.m.addConstr( self.g[s, s, v, 'node'] == 0 )
                for t in self.S:
                    self.m.addConstr( self.g[s, t, v, 'node'] + self.g[t, s, v, 'node'] <= 1 )

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

        if criteria_ordering[0] == 'C':
            for s in self.S:
                for t in self.S:
                    self.m.addConstr(self.z[s, t] <= 1)
            self.m.addConstr( self.z.sum() == len(self.S) - 1 )

        if criteria_ordering[0] == 'S':
            self.m.addConstr(self.s.sum() == 1)


    def _add_primary_location_constraint(self, primary_location):
        """Private function for internal use by the MACH2 algorithm."""
        self.m.addConstr( self.r[primary_location] == 1 )

    def _count_retrieved_solutions(self):
        """Private function for internal use by the MACH2 algorithm."""
        self.m.setParam(GRB.Param.SolutionNumber, 0)
        best_obj_val = self.m.PoolObjVal
        for e in range(self.m.SolCount):
            self.m.setParam(GRB.Param.SolutionNumber, e)
            if self.m.PoolObjVal - best_obj_val > 0.5:
                return e
        return self.m.SolCount

    def solve(self, logfile='', starting_nsols = 37, max_solutions=37888, threads=0, timelimit=0):
        """
        Solves the MACH2 optimization problem and retrieves the solutions.
        
        :param logfile: Path to the log file for Gurobi output.
        :param starting_nsols: Initial number of solutions to retrieve.
        :param max_solutions: Maximum number of solutions to retrieve.
        :param threads: Number of threads to use for optimization.
        :param timelimit: Amount of time the ILP will be run.
        :return: A SolutionSet containing the refined trees (Refinement).
        """
        while True:
            self.m.setParam(GRB.Param.MIPGap, 0)
            self.m.setParam(GRB.Param.PoolSearchMode, 2)
            self.m.setParam(GRB.Param.Threads, threads)
            self.m.setParam(GRB.Param.LogToConsole, 0)
            self.m.setParam(GRB.Param.LogFile, logfile)
            self.m.setParam(GRB.Param.PoolSolutions, starting_nsols)
            if timelimit > 0:
                self.m.setParam(GRB.Param.TimeLimit, timelimit)
            self.m.optimize()

            if self.m.SolCount == 0:
                raise ValueError("No feasible solution found within time limit (default: None).")
            elif self.m.status == GRB.TIME_LIMIT:
                if self.m.MIPGap == 0:
                    print("Time limit reached. Partial optimal solution space returned.")
                else:
                    print("Time limit reached. Solutions may be suboptimal.")
                break

            n = self._count_retrieved_solutions()
            if starting_nsols > n:
                self.successfully_run = True
                break
            elif starting_nsols * 2 > max_solutions:
                self.successfully_run = False
                print("Maximum number of solutions reached. Partial optimal solution space returned.")
                break
            else:
                starting_nsols *= 2
                self.m.reset(1)
        # print(n)
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
        
        solset = SolutionSet([Refinement.from_refinement_graph(self.tree, R) for R in Rs])
        if 'C' in self.criteria_ordering:
            inferred_comigs = sum([self.z[s,t].X for s in self.S for t in self.S])
            actual_comigs = [s.n_comigrations() for s in solset]
            if inferred_comigs == max(actual_comigs):
                return solset
            elif inferred_comigs == min(actual_comigs):
                return solset.filter(self.criteria_ordering)
            elif min(actual_comigs) - inferred_comigs < 2:
                return solset.filter(self.criteria_ordering)
            else:
                raise ValueError("Super exceptional case detected. Please open a GitHub issue in MACH2 repository and I'll implement this case.")
        else:
            return solset
