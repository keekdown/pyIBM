# file: $pyIBM/src/variable.py
# author: Olivier Mesnard (mesnardo@gwu.edu)
# BarbaGroup (lorenabarba.com)


import os

import numpy as np
import yaml

from case import Case
from parameters import Parameters
from mesh import Mesh
from matrix import Matrix


class BoundaryCondition:
	"""Defines a boundary condition by its type and values."""
	def __init__(self, value, N):
		"""Creates the boundary condition.
		
		Arguments
		---------
		value -- contains type and value at the boundary.
		N -- number of points on the boundary.
		"""
		self.bc_type = value[0]
		self.values = value[1] * np.ones(N, dtype=float)


class BoundaryConditions:
	"""Boundary conditions related to a variable."""
	def __init__(self, info_bc):
		"""Creates a list for each boundary with the boundary type and the boundary value array.
		
		Arguments
		---------
		info_bc -- information related to boundary conditions of a variable.
		"""
		for location, value in info_bc.iteritems():
			N = (Mesh.Nx if location in ['bottom', 'top'] else Mesh.Ny)
			setattr(self, location, BoundaryCondition(value, N))


class Variable:
	"""Creates a variable."""
	def __init__(self, name, info):
		"""Creates a variable parsing the file _infoFlow.yaml.

		Arguments
		---------
		name -- variable's name.
		info -- info about initial conditions and boundary conditions.
		skip_assemble -- boolean, if True, not assembling matrices related to the variable (default False).
		"""
		self.name = name
		# get initial conditions or read a data file
		self.get_initial_conditions(info['initialCondition'])
		# boundary conditions
		self.set_boundary_conditions(info['boundaryCondition'])

	def get_initial_conditions(self, info_ic):
		"""Gets the initial conditions.
		
		Arguments
		---------
		info_ic -- dictionary that contains the info related to the initial conditions.
		"""
		if Parameters.start == 0:
			self.field = ( info_ic
						 * np.ones(Mesh.Nx*Mesh.Ny, dtype=float) )
		else:
			self.read()

	def set_boundary_conditions(self, info_bc):
		"""Sets the boundary conditions.
		
		Arguments
		---------
		info_bc -- dictionary that contains the info related to the boundary conditions.
		"""
		self.bc = {}
		for location, value in info_bc.iteritems():
			N = (Mesh.Nx if location in ['bottom', 'top'] else Mesh.Ny)
			self.bc[location] = BoundaryCondition(value, N)

	def assemble_matrix(self, info):
		"""Assembles a matrix related to one variable.
		
		Arguments
		---------
		info -- info related to the matrix.
		"""
		setattr(self, info['type']+info['direction'],
				Matrix(self.bc, info['type'], info['scheme'],
					   info['direction']))

	def assemble_matrices(self):
		"""Assembles matrices related to one variable."""	
		with open(Case.path+'/_infoScheme.yaml', 'r') as infile:
			info_matrices = yaml.load(infile)
		for info_matrix in info_matrices[self.name]:
			if 'direction' not in info_matrix:
				info_matrix['direction'] = ''
			self.assemble_matrix(info_matrix)

	def write(self):
		"""Writes the variable field into a file."""
		with open(Case.path+'/'+str(Parameters.ite)+'/'+self.name+'.dat', 'w') as outfile:
			np.savetxt(outfile, np.c_[self.field], 
					   fmt='%.6f', delimiter='\t', 
					   header='%s - %d ites' % (self.name, Parameters.ite))
	
	def read(self):
		"""Reads the variable field from a file."""
		with open(Case.path+'/'+str(Parameters.ite)+'/'+self.name+'.dat', 'r') as infile:
			self.field = np.loadtxt(infile, 
									dtype=float, delimiter='\t')
