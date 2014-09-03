# file: $pyIBM/src/parameters.py
# author: Olivier Mesnard (mesnardo@gwu.edu)
# BarbaGroup (lorenabarba.com)


import yaml

from case import Case


class Parameters:
	"""Contains the parameters of the simulation."""
	def __init__(self):
		"""Parses the file _infoSolver.yaml and stores the parameters."""
		with open(Case.path+'/_infoSolver.yaml', 'r') as infile:
			info = yaml.load(infile)

		# parameters related to time
		Parameters.scheme = info['time']['scheme']
		Parameters.start = info['time']['start']
		Parameters.ite = Parameters.start
		Parameters.nt = info['time']['nt']
		Parameters.dt = info['time']['dt']
		Parameters.write_every = info['time']['writeEvery']
		
		# Reynolds number
		Parameters.Re = info['Re']
