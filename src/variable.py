# source: $pyIBM/src/variable.py
# Olivier Mesnard (mesnardo@gwu.edu)
# BarbaGroup (lorenabarba.com)

import os
import sys
import numpy as np
from scipy.sparse import *
import yaml

from case import Case
from mesh import Mesh
from solver import Solver
from matrix import *

class BoundaryConditions:
	'''Create the boundary condition with the type and value.
	Type: Dirichlet or Neummann.
	'''
	def __init__(self, info_bc):
		self.left = [info_bc['left'][0],info_bc['left'][1] * np.ones(Mesh.Ny, dtype=float)]
		self.right = [info_bc['right'][0],info_bc['right'][1] * np.ones(Mesh.Ny, dtype=float)]
		self.bottom = [info_bc['bottom'][0],info_bc['bottom'][1] * np.ones(Mesh.Nx, dtype=float)]
		self.top = [info_bc['top'][0],info_bc['top'][1] * np.ones(Mesh.Nx, dtype=float)]

class Variable:
	'''Create a variable.'''

	def __init__(self, name):
		self.name = name
		infile = open(Case.path+'/_infoFlow.yaml', 'r')
		info = yaml.load(infile)
		infile.close()
		if Solver.start == 0:
			self.field = info[self.name]['initialCondition']\
						* np.ones(Mesh.Nx*Mesh.Ny, dtype=float)
		else:
			self.field = np.empty(Mesh.Nx*Mesh.Ny, dtype=float)
			self.read()
		self.bc = BoundaryConditions(info[self.name]['boundaryCondition'])
	

	def assemble_matrix(self, mat_name, scheme='central', direction=None):
		'''Assemble a matrix related to a variable.'''

		setattr(self, mat_name, Matrix(self.bc, mat_name, scheme, direction))


	def read(self):
		'''Read the variable field from a file.'''

		infile = open(Case.path+'/'+str(Solver.ite)+'/'+self.name+'.dat', 'r')
		i = 0
		for line in infile:
			data = line.split()
			self.field[i] = float(data[0])
			i += 1
		infile.close()


	def write(self):
		'''Write the variable field into a file.'''

		if not os.path.isdir(Case.path+'/'+str(Solver.ite)):
			os.system('mkdir '+Case.path+'/'+str(Solver.ite))
		outfile = open(Case.path+'/'+str(Solver.ite)+'/'+self.name+'.dat', 'w')
		for i in xrange(len(self.field)):
			outfile.write(str(self.field[i])+'\n')
		outfile.close()
