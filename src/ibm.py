# file: $pyIBM/src/ibm.py
# author: Olivier Mesnard (mesnardo@gwu.edu)
# BarbaGroup (lorenabarba.com)


from math import sqrt

import numpy as np

from mesh import Mesh
from parameters import Parameters


def dh(r):
	"""Computes the delta function (Roma et al., 1999).
	
	Arguments
	---------
	r -- distance between Eulerian point and Lagrangian point.
	"""
	r = abs(r)
	if 0. <= r < 1.:
		return 0.125*(3-2*r+sqrt(1+4*r-4*r**2))
	elif 1. <= r < 2.:
		return 0.125*(5-2*r-sqrt(-7+12*r-4*r**2))
	else:
		return 0.0

def interpolation(u, body):
	"""Interpolates the velocity at the Lagrangian points.

	Arguments
	---------
	u -- component of the Eulerian velocity field.
	body -- Body object containing info about the immersed body.

	Returns
	-------
	uk -- component of the Lagangian velocity.
	"""
	uk = np.zeros(body.N, dtype=float)
	for k, (neighbor, x, y) in enumerate(zip(body.neighbor, body.x, body.y)):
		h = ( Mesh.x[neighbor%Mesh.Nx+1]
			- Mesh.x[neighbor%Mesh.Nx] )
		for i in xrange(-2,3):
			for j in xrange(-2,3):
				I = neighbor%Mesh.Nx + i
				J = neighbor/Mesh.Nx + j
				uk[k] += ( u[J*Mesh.Nx+I]
					     * dh((x - Mesh.x[I])/h)
					     * dh((y - Mesh.y[J])/h) )
	return uk

def distribution(fk, body):
	"""Spreads the Lagrangian force onto the Eulerian mesh grid.

	Arguments
	---------
	fk -- Lagrangian force.
	body -- Body object containing info about the immersed body.

	Returns
	-------
	f -- component of the Eulerian force field.
	"""
	f = np.zeros(Mesh.Nx*Mesh.Ny, dtype=float)

	for neighbor, x, y, f_k in zip(body.neighbor, body.x, body.y, fk):
		h = ( Mesh.x[neighbor%Mesh.Nx+1]
			- Mesh.x[neighbor%Mesh.Nx] )
		for i in xrange(-2,3):
			for j in xrange(-2,3):
				I = neighbor%Mesh.Nx + i
				J = neighbor/Mesh.Nx + j
				f[J*Mesh.Nx+I] += ( f_k
								  * dh((x - Mesh.x[I])/h)
								  * dh((y - Mesh.y[J])/h) )
	return f


def ibm(body, u, v):
	"""Immersed boundary method to update the Eulerian velocity field.
	
	Arguments
	---------
	body -- Body object containing info about the immersed body.
	u, v -- Eulerian velocity field.
	"""
	fx = np.zeros(Mesh.Nx*Mesh.Ny, dtype=float)
	fy = np.zeros(Mesh.Nx*Mesh.Ny, dtype=float)

	if body.is_moving:
		body.kinematics()
	
	dt = Parameters.dt
	N_ibm = 1
	body.cl,body.cd = 0.,0.
	for ite in xrange(N_ibm):
		body.u = interpolation(u, body)
		body.v = interpolation(v, body)
		body.fx[:] = (body.ud[:] - body.u[:]) / dt
		body.fy[:] = (body.vd[:] - body.v[:]) / dt

		for neighbor, fx_b, fy_b in zip(body.neighbor, body.fx, body.fy):
			body.cd += ( -2.*fx_b
					   * Mesh.dx[neighbor%Mesh.Nx]
					   * Mesh.dy[neighbor/Mesh.Nx] )
			body.cl += ( -2.*fy_b
					   * Mesh.dx[neighbor%Mesh.Nx]
					   * Mesh.dy[neighbor/Mesh.Nx] )

		fx += distribution(body.fx, body)
		fy += distribution(body.fy, body)

	u += dt * fx
	v += dt * fy
	
	return u, v
