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
	'''Creates the boundary condition with the type and value.
	Type of boundary condition: Dirichlet or Neummann.
	'''
	def __init__(self, info_bc):
		self.left = [info_bc['left'][0],info_bc['left'][1] * np.ones(Mesh.Ny, dtype=float)]
		self.right = [info_bc['right'][0],info_bc['right'][1] * np.ones(Mesh.Ny, dtype=float)]
		self.bottom = [info_bc['bottom'][0],info_bc['bottom'][1] * np.ones(Mesh.Nx, dtype=float)]
		self.top = [info_bc['top'][0],info_bc['top'][1] * np.ones(Mesh.Nx, dtype=float)]


class Variable:
	'''Creates a variable.'''
	def __init__(self, name, skip_assemble=False):
		self.name = name
		infile = open(Case.path+'/_infoFlow.yaml', 'r')
		info = yaml.load(infile)
		infile.close()	
		# get initial conditions or read a data file
		if Solver.start == 0:
			self.field = ( info[self.name]['initialCondition']
						 * np.ones(Mesh.Nx*Mesh.Ny, dtype=float) )
		else:
			self.read()
		# boundary conditions
		self.bc = BoundaryConditions(info[self.name]['boundaryCondition'])
		# assemble matrices
		if not skip_assemble:
			infile = open(Case.path+'/_infoScheme.yaml', 'r')
			info = yaml.load(infile)
			infile.close()
			for d in info[self.name]:
				if 'direction' not in d:
					d['direction'] = ''
				self.assemble_matrix(name=d['type'],
									 scheme=d['scheme'], 
									 direction=d['direction'])

	def assemble_matrix(self, name, scheme, direction):
		'''Assembles a matrix related to a variable,
		calling the class Matrix.
		'''
		setattr(self, name+direction, 
				Matrix(self.bc, name, scheme, direction))

	def write(self):
		'''Writes the variable field into a file.'''
		if not os.path.isdir(Case.path+'/'+str(Solver.ite)):
			os.system('mkdir '+Case.path+'/'+str(Solver.ite))
		with open(Case.path+'/'+str(Solver.ite)+'/'+self.name+'.dat', 'w') as file_name:
			np.savetxt(file_name, np.c_[self.field], 
					   fmt='%.6f', delimiter='\t', 
					   header='%s - %d ites' % (self.name, Solver.ite))
	
	def read(self):
		'''Reads the variable field from a file.'''
		with open(Case.path+'/'+str(Solver.ite)+'/'+self.name+'.dat', 'r') as file_name:
			self.field = np.loadtxt(file_name, 
									dtype=float, delimiter='\t')
