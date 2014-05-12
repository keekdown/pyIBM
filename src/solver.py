# source: $pyIBM/src/solver.py
# Olivier Mesnard (mesnardo@gwu.edu)
# BarbaGroup (lorenabarba.com)

import os
import numpy as np
import yaml

from mesh import Mesh
from case import Case

class Solver:
	'''Create the solver by reading the solver file
	in the case folder.
	'''
	def __init__(self):
		infile = open(Case.path+'/_infoSolver.yaml', 'r')
		info = yaml.load(infile)
		infile.close()
		Solver.scheme = info['time']['scheme']
		Solver.start = info['time']['start']
		Solver.ite = Solver.start
		Solver.nt = info['time']['nt']
		Solver.dt = info['time']['dt']
		Solver.write_every = info['time']['writeEvery']
		Solver.Re = info['Re']
