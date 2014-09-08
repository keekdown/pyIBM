#!/usr/env/bin python
# file: $pyIBM/pyIBM.py
# author: Olivier Mesnard (mesnardo@gwu.edu)
# BarbaGroup (lorenabarba.com)


import os
import sys

sys.path.insert(0,os.getcwd()+'/src')
from case import Case
from domain import Domain
from solver import Solver
from input_output import read_inputs


def main():
	"""Solves the Navier-Stokes equations in a two-dimensional domain
	with an immersed boundary method."""

	# parse the command-line
	args = read_inputs()

	# create the case
	Case(args.path)

	# create the computational domain
	Domain(is_show=args.mesh)
	
	# create the solver
	solver = Solver()

	# solve the Navier-Stokes equations
	solver.solve(body=Domain.body)
	

if __name__ == '__main__':
	print '\n\t----- pyIBM - START -----\n'
	main()
	print '\n\t----- pyIBM - END -----\n'
