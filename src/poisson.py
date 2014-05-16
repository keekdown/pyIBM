# source: $pyIBM/src/poisson.py
# Olivier Mesnard (mesnardo@gwu.edu)
# BarbaGroup (lorenabarba.com)

import numpy as np
import scipy.sparse.linalg as spla
import yaml
import pyamg

from case import Case


class Poisson:
	'''Creates a solver to solve a Poisson equation.'''
	solvers = {'cg':spla.cg, 'gmres':spla.gmres, 'bicg':spla.bicg}

	def __init__(self, var):
		'''Parses the file _infoSolver.'''
		infile = open(Case.path+'/_infoSolver.yaml', 'r')
		data = yaml.load(infile)
		infile.close()
		self.solver = Poisson.solvers[data['poisson']['solver']]
		self.tol = data['poisson']['tol']
		self.maxiter = data['poisson']['maxiter']
		self.M = None
		if data['poisson']['precond'] != None:
			self.ml = pyamg.smoothed_aggregation_solver(var.laplacian.mat)
			self.M = self.ml.aspreconditioner(cycle=data['poisson']['precond']['cycle'])
		self.iterations, self.residuals = [], []

	def solve(self, A, b, xi):
		'''Solves the Poisson equation,
		and stores the number of iterations and the residual.
		'''
		self.ite = 0
		x, info = self.solver(A, b, xi,
							  tol=self.tol, maxiter=self.maxiter, 
							  M=self.M, callback=self.iteration)
		self.iterations.append(self.ite)
		self.residuals.append(self.residual(A, x, b))
		return x

	def iteration(self, x):
		self.ite += 1
	
	def residual(self, A, x, b):
		'''Computes the residual using the L2-norm.'''
		if np.linalg.norm(b) != 0:
			return np.linalg.norm(A.dot(x)-b)/np.linalg.norm(b)
