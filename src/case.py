# source: $pyIBM/src/case.py
# Olivier Mesnard (mesnardo@gwu.edu)
# BarbaGroup (lorenabarba.com)
		
import os


class Case:
	"""Sets up information about the case."""
	def __init__(self,loc):
		"""Gets the path and name of the case.
		
		Arguments
		---------
		loc -- location of the case folder (root is the path of pyIBM.py)
		"""
		Case.path = os.getcwd()+'/'+loc
		Case.name = os.path.basename(os.path.normpath(Case.path))
