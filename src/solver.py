# source: $pyIBM/src/solver.py
# Olivier Mesnard (mesnardo@gwu.edu)
# BarbaGroup (lorenabarba.com)

import os
import numpy as np
import yaml

from case import Case


class Solver:
	"""Creates a solver."""
	def __init__(self):
		"""Parses the file _infoSolver.yaml to get solver info."""
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
