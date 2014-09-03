# file: $pyIBM/src/solver.py
# author: Olivier Mesnard (mesnardo@gwu.edu)
# BarbaGroup (lorenabarba.com)


import os

import numpy as np
import yaml

from case import Case
from parameters import Parameters
from mesh import Mesh
from variable import Variable
from poisson import Poisson
from ibm import ibm
from operations import grad, lap
import time_info


class Solver:
	"""Creates a solver."""
	def __init__(self):
		"""Initializes the solver.
		Creates the flow variables, 
		assembles the discretized operator matrices
		and initializes the Poisson solver.
		"""
		# read the simulation parameters
		params = Parameters()

		print '\n{Assembling matrices}'
		tic = time_info.start()
		Solver.u = Variable('u')
		Solver.v = Variable('v')
		Solver.p = Variable('p')
		time_info.stop(tic, 'Assembling matrices')

		Solver.poisson = Poisson(Solver.p)

	def solve(self, body=None):
		# copy the variables for readability
		u = Solver.u
		v = Solver.v
		p = Solver.p

		# copy some parameters for readability
		ite = Parameters.ite
		start = Parameters.start
		nt = Parameters.nt
		dt = Parameters.dt
		Re = Parameters.Re

		# open the file that contains force coefficients
		outfile = open(Case.path+'/forceCoeffs.dat',
					   ('w' if Parameters.start == 0 else 'a'))

		# time-integration using Euler method
		tic = time_info.start()
		while ite < start + nt:
			ite += 1
			print '\nIteration %d - Time = %f' % (ite, ite*dt)
			# updates velocity field without the pressure
			u.field += dt * (1./Re*lap(u) - u.field*grad(u, 'x') - v.field*grad(u, 'y'))
			v.field += dt * (1./Re*lap(v) - u.field*grad(v, 'x') - v.field*grad(v, 'y'))

			# immersed boundary method
			if body:
				ibm(body, u, v)
				outfile.write('%f \t %f \t %f' % (ite*dt, body.cl, body.cd))

			# solve Poisson equation for pressure
			b = 1./dt * (grad(u, 'x') + grad(v, 'y'))
			p.field = Solver.poisson.solve(p.laplacian.mat, b-p.laplacian.bc_vect, p.field)

			print '{Poisson} Number of iterations: %d' % Solver.poisson.ite

			# update velocity field
			u.field -= dt * grad(p, 'x')
			v.field -= dt * grad(p, 'y')

			# write variable fields
			if ite%Parameters.write_every == 0:
				print '\n{Writting results}'
				u.write()
				v.write()
				p.write()
		
		time_info.stop(tic, 'DONE')
		
		# close the file containing force coefficients
		outfile.close()
