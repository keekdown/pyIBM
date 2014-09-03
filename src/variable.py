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
	def __init__(self, name, skip_assemble=False):
		"""Creates a variable parsing the file _infoFlow.yaml.

		Arguments
		---------
		name -- variable's name.
		skip_assemble -- boolean, if True, not assembling matrices related to the variable (default False).
		"""
		self.name = name
	
		# parse flow file with yaml
		with open(Case.path+'/_infoFlow.yaml', 'r') as infile:
			info = yaml.load(infile)
		
		# get initial conditions or read a data file
		if Parameters.start == 0:
			self.field = ( info[self.name]['initialCondition']
						 * np.ones(Mesh.Nx*Mesh.Ny, dtype=float) )
		else:
			self.read()

		# boundary conditions
		self.set_boundary_conditions(info[self.name]['boundaryCondition'])
		#self.bc = BoundaryConditions(info[self.name]['boundaryCondition'])
		
		# assemble matrices
		if not skip_assemble:
			# parse schem file using yaml
			with open(Case.path+'/_infoScheme.yaml', 'r') as infile:
				info = yaml.load(infile)
			for d in info[self.name]:
				self.assemble_matrix(name=d['type'],
									 scheme=d['scheme'], 
									 direction=(d['direction'] if 'direction' in d else ''))

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

	def assemble_matrix(self, name, scheme, direction):
		"""Assembles a matrix related to a variable, calling the class Matrix.
		
		Arguments
		---------
		name -- variable's name.
		scheme -- numerical scheme for discretization.
		direction -- gradient's direction.
		"""
		setattr(self, name+direction, 
				Matrix(self.bc, name, scheme, direction))

	def write(self):
		"""Writes the variable field into a file."""
		if not os.path.isdir(Case.path+'/'+str(Parameters.ite)):
			os.system('mkdir '+Case.path+'/'+str(Parameters.ite))
		with open(Case.path+'/'+str(Parameters.ite)+'/'+self.name+'.dat', 'w') as file_name:
			np.savetxt(file_name, np.c_[self.field], 
					   fmt='%.6f', delimiter='\t', 
					   header='%s - %d ites' % (self.name, Parameters.ite))
	
	def read(self):
		"""Reads the variable field from a file."""
		with open(Case.path+'/'+str(Parameters.ite)+'/'+self.name+'.dat', 'r') as file_name:
			self.field = np.loadtxt(file_name, 
									dtype=float, delimiter='\t')
