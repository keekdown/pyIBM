# file: $pyIBM/src/case.py
# author: Olivier Mesnard (mesnardo@gwu.edu)
# BarbaGroup (lorenabarba.com)
		

import os

class Case:
	"""Sets up information about the case."""
	def __init__(self, path):
		"""Gets the path and name of the case.
		
		Arguments
		---------
		path -- location of the case folder (root is the path of pyIBM.py).
		"""
		Case.path = os.getcwd() + '/' + path
		Case.name = os.path.basename(os.path.normpath(Case.path))
		Case.images  = Case.path + '/images'

		# create sub-folder containing images
		if not os.path.isdir(Case.images):
			os.system('mkdir ' + Case.images)
