# script: $pyIBM/scripts/clean.py
# Olivier Mesnard (mesnardo@gwu.edu)
# BarbaGroup (lorenabarba.com)

import os
import sys


def main(arg):
	"""Script to clean a given case."""
	pwd = os.getcwd()
	case_name = arg[1]
	case_path = pwd+'/'+case_name+'/'

	print '[info]: cleaning solution folders in ', case_name
	for i in xrange(10):
		os.system('rm -rf '+case_path+str(i)+'*')

	print '[info]: cleaning images in ', case_name
	os.system('rm -rf '+case_path+'images')
	os.system('rm -f '+case_path+'*.dat')
	os.system('rm -f '+case_path+'*.png')

if __name__ == '__main__':
	print '\n\t----- START - CLEANING -----\n'
	main(sys.argv)
	print '\n\t----- END - CLEANING -----\n'
