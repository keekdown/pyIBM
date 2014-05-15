# script: $pyIBM/scripts/geometry_cylinder.py
# Olivier Mesnard (mesnardo@gwu.edu)
# BarbaGroup (lorenabarba.com)

import os
import sys
import argparse
from math import *
import numpy as np

def main():
	'''Creates a cylinder with the same spatial distribution than
	the uniform where it will be immersed.
	'''
	
	# list of command-line arguments
	parser = argparse.ArgumentParser(description='Creates a cylinder with the same spatial 
									 distribution than the uniform grid where it will be immersed')
	parser.add_argument('-p', '--path', dest='path', 
						help='path of the case folder', type=str)
	parser.add_argument('-b', '--boundary', dest='boundary', 
						help='boundaries of the uniform grid', 
						nargs='+', type=float)
	parser.add_argument('-n', '--n', dest='n', 
						help='number of points in the x and y directions in the uniform region', 
						nargs='+', type=float)
	parser.add_argument('-c', '--center', dest='center', 
						help='center of the cylinder', 
						nargs='+', type=float, default=[0, 0])
	parser.add_argument('-r', '--radius', dest='radius', 
						help='radius of the cylinder', 
						type=float, default=1.0)
	args = parser.parse_args()

	# gets the path of the case
	pwd = os.getcwd()
	case_name = args.path
	case_path = pwd+'/'+case_name

	# gets boundaries and number of points
	xmin, xmax = args.boundary[0], args.boundary[1]
	Nx, Ny = args.n[0], args.n[1]

	# gets radius and center of the cylinder
	R = args.radius
	x_center, y_center = args.center[0], args.center[1]

	# computes coordinates of the cylinder
	N = int(2*pi*R*(Nx-1)/(xmax-xmin))
	x = np.empty(N, dtype=float)
	y = np.empty(N, dtype=float)
	for k in xrange(N):
		x[k] = x_center + R*cos(2*pi*k/N)
		y[k] = y_center + R*sin(2*pi*k/N)

	# writes coordinates in a file in the case folder
	outfile = open(case_path+'/cylinder.bdy', 'w')
	outfile.write(str(N)+'\n')
	for k in xrange(N):
		outfile.write(str(x[k])+'\t'+str(y[k])+'\n')
	outfile.close()

if __name__ == '__main__':
	print '\n\t----- pyIBM - START - cylinder geometry -----\n'
	main()
	print '\n\t----- pyIBM - END - cylinder geometry -----\n'
