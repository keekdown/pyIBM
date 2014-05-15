# source: $pyIBM/src/poisson.py
# Olivier Mesnard (mesnardo@gwu.edu)
# BarbaGroup (lorenabarba.com)

import numpy as np
import scipy.sparse.linalg as spla
import yaml
import pyamg

from case import Case

class Poisson:
	'''Create Poisson solvers.'''

	solvers = {'cg':spla.cg, 'gmres':spla.gmres, 'bicg':spla.bicg}

	def __init__(self, var):
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
		'''Solve Poisson equation,
		and store the number of iterations and residual.
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
		if np.linalg.norm(b) != 0:
			return np.linalg.norm(A.dot(x)-b)/np.linalg.norm(b)
