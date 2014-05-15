# source: $pyIBM/src/matrix.py
# Olivier Mesnard (mesnardo@gwu.edu)
# BarbaGroup (lorenabarba.com)

import os
import sys
import numpy as np
from scipy.sparse import *
import yaml

from mesh import Mesh


class Matrix:
	'''Create a matrix (laplacian or gradient), given a scheme.
	The direction (x or y) to compute the gradient matrix.
	'''
	def __init__(self, bc, name, scheme, direction):
		self.bc = bc
		self.mat, self.bc_vect = getattr(self, 'set_'+name)(direction, scheme)
	
	def set_gradient(self, direction, scheme='central'):
		'''Create gradient matrix,
		given a scheme and a direction.
		'''
		Nx, Ny = Mesh.Nx, Mesh.Ny
		A = []
		B = np.zeros(Nx*Ny, dtype=float)
		if direction == 'x':
			N = Nx
			d = Mesh.dx
			bc1, bc2 = self.bc.left, self.bc.right
			p = 1
		elif direction == 'y':
			N = Ny
			d = Mesh.dy
			bc1, bc2 = self.bc.bottom, self.bc.top
			p = Nx
		if scheme == 'backward':
			for i in xrange(Nx*Ny):
				if direction == 'x':
					I, J = i%Nx, i/Nx
				elif direction =='y':
					I, J = i/Nx, i%Nx
				#bc1
				if I == 0:
					if bc1[0] == 'dirichlet':
						A.append([i, i, 1./d[I]])
						B[i] += -bc1[1][J]/d[I]
					elif bc1[0] == 'neumann':
						B[i] += bc1[1][J]
				else:
					A.append([i, i-p, -1./d[I-1]])
					A.append([i, i, 1./d[I-1]])
		elif scheme == 'central':
			for i in xrange(Nx*Ny):
				if direction == 'x':
					I, J = i%Nx, i/Nx
				elif direction =='y':
					I, J = i/Nx, i%Nx
				#bc1
				if I == 0:
					if bc1[0] == 'dirichlet':
						A.append([i, i+p, 0.5/d[I]])
						B[i] += -0.5*bc1[1][J]/d[I]
					elif bc1[0] == 'neumann':
						B[i] += bc1[1][J]
				#bc2
				elif I == N-1:
					if bc2[0] == 'dirichlet':
						A.append([i, i-p, -0.5/d[I-1]])
						A.append([i, i, 0.5*(1./d[I-1]-1./d[I])])
						B[i] += 0.5*bc2[1][J]/d[I]
					elif bc2[0] == 'neumann':
						B[i] += bc2[1][J]
				else:
					A.append([i, i-p, -0.5/d[I-1]])
					A.append([i, i, 0.5*(1./d[I-1]-1./d[I])])
					A.append([i, i+p, 0.5/d[I]])
		val, row, col = [x[2] for x in A], [x[0] for x in A], [x[1] for x in A]
		
		return csr_matrix((val, (row, col)), shape=(Nx*Ny, Nx*Ny), dtype=float), B

	def set_laplacian(self, direction=None, scheme='central'):
		'''Create Laplacian matrix,
		using central difference scheme.
		'''
		Nx, Ny = Mesh.Nx, Mesh.Ny
		dx, dy = Mesh.dx, Mesh.dy
		bc = self.bc
		A = []
		B = np.zeros(Nx*Ny, dtype=float)
		for i in xrange(Nx*Ny):
			I, J = i%Nx, i/Nx
			#left
			if I == 0:
				A.append([i, i+1, 1./dx[I]**2])
				A.append([i, i, -2./dx[I]**2])
				if bc.left[0] == 'dirichlet':
					B[i] += bc.left[1][J]/dx[I]**2
				elif bc.left[0] == 'neumann':
					A.append([i, i, 1./dx[I]**2])
					B[i] += -bc.left[1][J]/dx[I]
			#right
			if I == Nx-1:
				A.append([i, i-1, 1./dx[I]**2])
				A.append([i, i, -2./dx[I]**2])
				if bc.right[0] == 'dirichlet':
					B[i] += bc.right[1][J]/dx[I]**2
				elif bc.right[0] == 'neumann':
					A.append([i, i, 1./dx[I]**2])
					B[i] += bc.right[1][J]/dx[I]
			#bottom
			if i < Nx:
				A.append([i, i+Nx, 1./dy[J]**2])
				A.append([i, i, -2./dy[J]**2])
				if bc.bottom[0] == 'dirichlet':
					B[i] += bc.bottom[1][I]/dy[J]**2
				elif bc.bottom[0] == 'neumann':
					A.append([i, i, 1./dy[J]**2])
					B[i] += -bc.bottom[1][I]/dy[J]
			#top
			if i >= Nx*(Ny-1):
				A.append([i, i-Nx, 1./dy[J]**2])
				A.append([i, i, -2./dy[J]**2])
				if bc.top[0] == 'dirichlet':
					B[i] += bc.top[1][I]/dy[J]**2
				elif bc.top[0] == 'neumann':
					A.append([i, i, 1./dy[J]**2])
					B[i] += bc.top[1][I]/dy[J]
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
