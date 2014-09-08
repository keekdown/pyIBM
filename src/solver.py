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
from input_output import timer_start, timer_stop


class Solver:
	"""Creates a solver."""
	def __init__(self, skip_assemble=False, skip_poisson=False):
		"""Initializes the solver.
		
		Creates the flow variables, 
		assembles the discretized operator matrices
		and initializes the Poisson solver.
		"""
		# read the simulation parameters
		Parameters()

		# initialize flow variables
		self.initialize_variables()

		print '\n{Assembling matrices}'
		tic = timer_start()
		if not skip_assemble:
			self.assemble_matrices()
		else:
			print '--> skipped'
		timer_stop(tic, info='Assembling matrices')

		if not skip_poisson:
			Solver.poisson = Poisson(Solver.p)
		else:
			print '--> Poisson skipped'

	def initialize_variables(self):
		with open(Case.path+'/_infoFlow.yaml', 'r') as infile:
			info_variables = yaml.load(infile)
		for name, info in info_variables.iteritems():
			setattr(Solver, name, Variable(name, info))
	
	def assemble_matrices(self):
		Solver.u.assemble_matrices()
		Solver.v.assemble_matrices()
		Solver.p.assemble_matrices()

	def solve(self, body=None):

		# open the file that contains force coefficients
		outfile = open(Case.path+'/forceCoeffs.dat',
					   ('w' if Parameters.start == 0 else 'a'))

		# time-integration using Euler method
		tic = timer_start()
		while Parameters.ite < Parameters.start + Parameters.nt:
			Parameters.ite += 1
			print '\nIteration %d - Time = %f' \
				  % (Parameters.ite, Parameters.ite*Parameters.dt)

			# get intermediate velocity
			self.intermediate_velocity()

			# immersed boundary method
			self.immersed_boundary_method(body)

			# solve Poisson equation for pressure
			self.solve_poisson()

			# update velocity field
			self.update_velocity()

			# write variable fields
			if Parameters.ite%Parameters.write_every == 0:
				print '\n{Writting results}'
				if not os.path.isdir(Case.path+'/'+str(Parameters.ite)):
					os.system('mkdir '+Case.path+'/'+str(Parameters.ite))
				Solver.u.write()
				Solver.v.write()
				Solver.p.write()
			# write force coefficients
			outfile.write('%f \t %f \t %f' 
						  % (Parameters.ite*Parameters.dt, body.cl, body.cd))
		
		timer_stop(tic, info='DONE')
		# close the file containing force coefficients
		outfile.close()

	def intermediate_velocity(self):
		"""Computes the intermediate velocity field.
		
		Solves the Navier-Stokes equations,
		without the pressure gradient and boundardy forces,
		using Euler's method.
		"""
		Solver.u.field += Parameters.dt * (
							1./Parameters.Re*lap(Solver.u)
							- Solver.u.field*grad(Solver.u, 'x')
							- Solver.v.field*grad(Solver.u, 'y'))
		Solver.v.field += Parameters.dt * (
							1./Parameters.Re*lap(Solver.v)
							- Solver.u.field*grad(Solver.v, 'x')
							- Solver.v.field*grad(Solver.v, 'y'))

	def immersed_boundary_method(self, body=None):
		"""Udpates the velocity field 
		with the effect of the immersed boundary.
		
		Arguments
		---------
		body -- immersed boundary (default None).
		"""
		if body:
			Solver.u.field, Solver.v.field = ibm(body, Solver.u.field, Solver.v.field)
	
	def solve_poisson(self):
		"""Solves the Poisson equation for the pressure field."""
		b = 1./Parameters.dt * (grad(Solver.u, 'x') + grad(Solver.v, 'y'))
		Solver.p.field = Solver.poisson.solve(Solver.p.laplacian.mat, b-Solver.p.laplacian.bc_vect, Solver.p.field)
		print '{Poisson} Number of iterations: %d' \
			  % Solver.poisson.ite

	def update_velocity(self):
		"""Updates the velocity field with pressure correction
		to make the field divergence-free.
		"""
		Solver.u.field -= Parameters.dt * grad(Solver.p, 'x')
		Solver.v.field -= Parameters.dt * grad(Solver.p, 'y')
