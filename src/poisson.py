# file: $pyIBM/src/poisson.py
# author: Olivier Mesnard (mesnardo@gwu.edu)
# BarbaGroup (lorenabarba.com)


import numpy as np
import scipy.sparse.linalg as spla
import yaml


from case import Case

class Poisson:
    """Solver to solve Poisson equation."""
    solvers = {'cg': spla.cg, 'gmres': spla.gmres, 'bicg': spla.bicg, 'minres': spla.minres}

    def __init__(self, p):
        """Parses the file _infoSolver to initialize the solver.
        
        Arguments
        ---------
        p -- field variable (Variable object) of the Poisson equation.
        """
        with open(Case.path+'/_infoSolver.yaml', 'r') as infile:
            info = yaml.load(infile)['poisson']
        self.solver = Poisson.solvers[info['solver']]
        self.tol = info['tol']
        self.maxiter = info['maxiter']
        self.M = None
        if 'precond' in info:
              import pyamg
              ml = pyamg.smoothed_aggregation_solver(p.laplacian.mat)
              self.M = ml.aspreconditioner(cycle=info['precond']['cycle'])
        self.iterations, self.residuals = [], []

    def solve(self, A, b, xi):
        """Solves the Poisson equation.

        Arguments
        ---------
        A -- discretized Poisson's matrix.
        b -- right hand-side of the discretized Poisson's equation.
        xi -- initial solution.

        Returns
        -------
        x -- solution of the Poisson's equation.
        """
        self.ite = 0
        x, info = self.solver(A, b, xi,
                              tol=self.tol, maxiter=self.maxiter, 
                              M=self.M, callback=self.iteration)
        self.iterations.append(self.ite)
        self.residuals.append(self.residual(A, x, b))
        return x

    def iteration(self, x):
        """Increments the number of iteration by 1.
        
        Arguments
        ---------
        x -- solution after the iteration (mandatory as argument but useless).
        """
        self.ite += 1
    
    def residual(self, A, x, b):
        """Returns the residual using the L2-norm.
        
        Arguments
        ---------
        A -- discretized Poisson's matrix.
        x -- solution of Poisson's equation.
        b -- right hand-side of Poisson's equation.
        """
        if np.linalg.norm(b) != 0:
            return np.linalg.norm(A.dot(x)-b)/np.linalg.norm(b)
