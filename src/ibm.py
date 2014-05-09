# source: $pyIBM/src/ibm.py
# Olivier Mesnard (mesnardo@gwu.edu)
# BarbaGroup (lorenabarba.com)

from math import *
import numpy as np

from mesh import Mesh

def dh(r):
	'''Compute the delta function,
	using the formulation from Roma et al. (1999).
	'''
	r = abs(r)
	if (0 <= r < 1):
		return 0.125*(3-2*r+sqrt(1+4*r-4*r**2))
	elif (1 <= r < 2):
		return 0.125*(5-2*r-sqrt(-7+12*r-4*r**2))
	else:
		return 0.0

def interpolation(u,body):
	'''Interpolate the velocity at the Lagrangian points,
	using values of the Eulerian of the neighbors.
	'''
	uk = np.zeros(body.N,dtype=float)
	for k in range(body.N):
		h = Mesh.x[body.neighbor[k]%Mesh.Nx+1]\
			- Mesh.x[body.neighbor[k]%Mesh.Nx]
		for i in range(-2,3):
			for j in range(-2,3):
				I = body.neighbor[k]%Mesh.Nx+i
				J = body.neighbor[k]/Mesh.Nx+j
				uk[k] += u[J*Mesh.Nx+I]\
						*dh((body.x[k]-Mesh.x[I])/h)\
						*dh((body.y[k]-Mesh.y[J])/h)
	return uk

def distribution(fk,body):
	'''Spread the Lagrangian force
	onto the Eulerian grid.
	'''
	f = np.zeros(Mesh.Nx*Mesh.Ny,dtype=float)
	for k in range(body.N):
		h = Mesh.x[body.neighbor[k]%Mesh.Nx+1]\
			- Mesh.x[body.neighbor[k]%Mesh.Nx]
		for i in range(-2,3):
			for j in range(-2,3):
				I = body.neighbor[k]%Mesh.Nx+i
				J = body.neighbor[k]/Mesh.Nx+j
				f[J*Mesh.Nx+I] += fk[k]\
								*dh((body.x[k]-Mesh.x[I])/h)\
								*dh((body.y[k]-Mesh.y[J])/h)
	return f
