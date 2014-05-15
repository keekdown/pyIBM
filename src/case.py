# source: $pyIBM/src/case.py
# Olivier Mesnard (mesnardo@gwu.edu)
# BarbaGroup (lorenabarba.com)
		
import os


class Case:
	'''Define the name and the path of the case.'''
	def __init__(self,loc):
		Case.path = os.getcwd()+'/'+loc
		Case.name = os.path.basename(os.path.normpath(Case.path))
