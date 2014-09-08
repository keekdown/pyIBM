# file: $pyIBM/scripts/clean_case.py
# author: Olivier Mesnard (mesnardo@gwu.edu)
# BarbaGroup (lorenabarba.com)


import os
import sys

sys.path.insert(0,'./src')
from case import Case
from input_output import read_inputs


def main():
	"""Script to clean a given case."""
	args = read_inputs()
	Case(args.path)

	print 'case path: %s' % Case.path

	print '[info]: cleaning mesh and force coefficients'
	os.system('rm -f '+Case.path+'/*.dat')

	print '[info]: cleaning solution folders'
	for i in xrange(10):
		os.system('rm -rf '+Case.path+'/'+str(i)+'*')

	print '[info]: cleaning images'
	os.system('rm -rf '+Case.images)

if __name__ == '__main__':
	print '\n\t----- START - CLEANING -----\n'
	main()
	print '\n\t----- END - CLEANING -----\n'
