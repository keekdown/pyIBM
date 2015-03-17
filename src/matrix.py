# file: $pyIBM/src/matrix.py
# author: Olivier Mesnard (mesnardo@gwu.edu)
# BarbaGroup (lorenabarba.com)

import os
import sys

import numpy as np
from scipy.sparse import *
import yaml

from mesh import Mesh


class Matrix:
	"""Class to define a matrix."""
	def __init__(self, bc, name, scheme, direction):
		"""Creates a Laplacian or gradient matrix associated to a variable.

		Arguments
		---------
		bc -- boundary conditions.
		name -- name of the matrix.
		scheme -- numerical scheme to differentiate the gradient or Laplacian.
		direction -- gradient's direction.
		"""
		self.bc = bc
		# call the right function with respect to the name
		self.mat, self.bc_vect = getattr(self, 'set_'+name)(direction, scheme)
	
	def set_gradient(self, direction, scheme='central'):
		"""Creates a gradient matrix and its boundary condition vector.
		
		Arguments
		---------
		direction -- gradient's direction.
		scheme -- numerical scheme for discretization (default 'central')
		"""
		Nx, Ny = Mesh.Nx, Mesh.Ny
		A = []
		B = np.zeros(Nx*Ny, dtype=float)

		dct = {'x': {'N': Mesh.Nx, 'delta': Mesh.dx,
					 'bc1': self.bc['left'], 'bc2': self.bc['right'], 'p': 1},
			   'y': {'N': Mesh.Ny, 'delta': Mesh.dy,
			   		 'bc1': self.bc['bottom'], 'bc2': self.bc['top'], 'p': Mesh.Nx}}

		N = dct[direction]['N']
		d = dct[direction]['delta']
		bc1, bc2 = dct[direction]['bc1'], dct[direction]['bc2']
		p = dct[direction]['p']

		if scheme == 'backward':
			for i in range(Nx*Ny):
				if direction == 'x':
					I, J = i%Nx, i/Nx
				elif direction =='y':
					I, J = i/Nx, i%Nx
				#bc1
				if I == 0:
					if bc1.bc_type == 'dirichlet':
						A.append([i, i, 1./d[I]])
						B[i] += -bc1.values[J]/d[I]
					elif bc1.bc_type == 'neumann':
						B[i] += bc1.values[J]
				else:
					A.append([i, i-p, -1./d[I-1]])
					A.append([i, i, 1./d[I-1]])
		elif scheme == 'central':
			for i in range(Nx*Ny):
				if direction == 'x':
					I, J = i%Nx, i/Nx
				elif direction =='y':
					I, J = i/Nx, i%Nx
				#bc1
				if I == 0:
					if bc1.bc_type == 'dirichlet':
						A.append([i, i+p, 0.5/d[I]])
						B[i] += -0.5*bc1.values[J]/d[I]
					elif bc1.bc_type == 'neumann':
						B[i] += bc1.values[J]
				#bc2
				elif I == N-1:
					if bc2.bc_type == 'dirichlet':
						A.append([i, i-p, -0.5/d[I-1]])
						A.append([i, i, 0.5*(1./d[I-1]-1./d[I])])
						B[i] += 0.5*bc2.values[J]/d[I]
					elif bc2.bc_type == 'neumann':
						B[i] += bc2.values[J]
				else:
					A.append([i, i-p, -0.5/d[I-1]])
					A.append([i, i, 0.5*(1./d[I-1]-1./d[I])])
					A.append([i, i+p, 0.5/d[I]])
		
		val, row, col = [x[2] for x in A], [x[0] for x in A], [x[1] for x in A]
		print(len(val),len(row),len(col),Nx*Ny)
		return csr_matrix((val, (row, col)), shape=(Nx*Ny, Nx*Ny), dtype=float), B

	def set_laplacian(self, direction=None, scheme='central'):
		"""Creates a Laplacian matrix and its boundary condition vector.
		
		Arguments
		---------
		direction -- None.
		scheme -- numerical scheme for discretization (default 'central')
		"""
		Nx, Ny = Mesh.Nx, Mesh.Ny
		dx, dy = Mesh.dx, Mesh.dy
		bc = self.bc
		A = []
		B = np.zeros(Nx*Ny, dtype=float)
		for i in range(Nx*Ny):
			I, J = i%Nx, i/Nx
			#left
			if I == 0:
				A.append([i, i+1, 1./dx[I]**2])
				A.append([i, i, -2./dx[I]**2])
				if bc['left'].bc_type == 'dirichlet':
					B[i] += bc['left'].values[J]/dx[I]**2
				elif bc['left'].bc_type == 'neumann':
					A.append([i, i, 1./dx[I]**2])
					B[i] += -bc['left'].values[J]/dx[I]
			#right
			if I == Nx-1:
				A.append([i, i-1, 1./dx[I]**2])
				A.append([i, i, -2./dx[I]**2])
				if bc['right'].bc_type == 'dirichlet':
					B[i] += bc['right'].values[J]/dx[I]**2
				elif bc['right'].bc_type == 'neumann':
					A.append([i, i, 1./dx[I]**2])
					B[i] += bc['right'].values[J]/dx[I]
			#bottom
			if i < Nx:
				A.append([i, i+Nx, 1./dy[J]**2])
				A.append([i, i, -2./dy[J]**2])
				if bc['bottom'].bc_type == 'dirichlet':
					B[i] += bc['bottom'].values[I]/dy[J]**2
				elif bc['bottom'].bc_type == 'neumann':
					A.append([i, i, 1./dy[J]**2])
					B[i] += -bc['bottom'].values[I]/dy[J]
			#top
			if i >= Nx*(Ny-1):
				A.append([i, i-Nx, 1./dy[J]**2])
				A.append([i, i, -2./dy[J]**2])
				if bc['top'].bc_type == 'dirichlet':
					B[i] += bc['top'].values[I]/dy[J]**2
				elif bc['top'].bc_type == 'neumann':
					A.append([i, i, 1./dy[J]**2])
					B[i] += bc['top'].values[I]/dy[J]
			#point not on left and not on right
			if I != 0 and I != Nx-1:
				A.append([i, i-1, 2./dx[I-1]/(dx[I-1]+dx[I])])
				A.append([i, i, -2./dx[I-1]/dx[I]])
				A.append([i, i+1, 2./dx[I]/(dx[I-1]+dx[I])])
			#point not on bottom and not on top
			if i >= Nx and i < Nx*(Ny-1):
				A.append([i, i-Nx, 2./dy[J-1]/(dy[J-1]+dy[J])])
				A.append([i, i, -2./dy[J-1]/dy[J]])
				A.append([i, i+Nx, 2./dy[J]/(dy[J-1]+dy[J])])
		
		val, row, col = [x[2] for x in A], [x[0] for x in A], [x[1] for x in A]
		
		return csr_matrix((val, (row, col)), shape=(Nx*Ny, Nx*Ny), dtype=float), B
