# script: $pyIBM/scripts/geometry_cylinder.py
# Olivier Mesnard (mesnardo@gwu.edu)
# BarbaGroup (lorenabarba.com)

import os
import sys
from math import *
import numpy as np

def main(arg):

	pwd = os.getcwd()
	case_name = arg[1]
	case_path = pwd+'/'+case_name

	xmin, xmax = -5.0, 10.0
	Nx, Ny = 120, 120

	R = 0.5
	x_center, y_center = 0.0, 0.0
	N = int(2*pi*R*(Nx-1)/(xmax-xmin))
	x = np.empty(N, dtype=float)
	y = np.empty(N, dtype=float)
	for k in xrange(N):
		x[k] = x_center + R*cos(2*pi*k/N)
		y[k] = y_center + R*sin(2*pi*k/N)
	outfile = open(case_path+'/cylinder.bdy', 'w')
	outfile.write(str(N)+'\n')
	for k in xrange(N):
		outfile.write(str(x[k])+'\t'+str(y[k])+'\n')
	outfile.close()

if __name__ == '__main__':
	print '\n\t----- pyIBM - START - cylinder geometry -----\n'
	main(sys.argv)
	print '\n\t----- pyIBM - END - cylinder geometry -----\n'
