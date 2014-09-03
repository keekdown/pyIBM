#!/usr/env/bin python

##############################################
# pyIBM - Immersed Boundary Method in Python #
# Olivier Mesnard (mesnardo@gwu.edu)         #
# BarbaGroup (lorenabarba.com)               #
##############################################

import os
import sys
import numpy as np

sys.path.insert(0,os.getcwd()+'/src')
from case import Case
from mesh import Mesh
from body import Body
from solver import Solver
from variable import Variable
from matrix import Matrix
from operations import grad, lap
from ibm import ibm
from poisson import Poisson

import time_info

def main(arg):

	Case(arg[0])

	# generate the mesh
	mesh = Mesh()

	# generate the immersed boundary
	if Mesh.is_body:
		body = Body()

	if '--mesh' in arg:
		mesh.plot(body if Mesh.is_body else None, is_show=True)
		sys.exit(0)
	else:
		mesh.plot(body if Mesh.is_body else None, is_show=False)

	# create the solver
	solver = Solver()

	# solve the Navier-Stokes equations
	solver.solve()
	

if __name__ == '__main__':
	print '\n\t----- pyIBM - START -----\n'
	main(sys.argv[1:])
	print '\n\t----- pyIBM - END -----\n'
