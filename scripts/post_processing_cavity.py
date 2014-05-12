# script: $pyIBM/scripts/post_processing_cavity.py
# Olivier Mesnard (mesnardo@gwu.edu)
# BarbaGroup (lorenabarba.com)

import os
import sys
import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0,'./src')
from case import *
from mesh import *
from solver import *
from variable import *

def main(arg):
	'''Plot velocity at center lines,
	with comparison Ghia et al. (1982) experimental data.
	'''
	Case(arg[1])

	mesh = Mesh()

	Solver()

	u = Variable('u')
	v = Variable('v')
	p = Variable('p')

	print '{Post-processing}: Comparison with Ghia et al. (1982)'
	Solver.ite = Solver.start + Solver.nt
	u.read()
	v.read()
	p.read()
	
	infile = open(os.getcwd()+'/resources/ghia_et_al_1982/u_vertical_line.dat', 'r')
	y_vl_ghia, u_vl_ghia = [], []
	index, row = 0, 0
	for line in infile:
		data = line.split()
		if index==0:
			for i in xrange(1, len(data)):
				if float(data[i]) == Solver.Re:
					row = i
				index += 1
		else:
			y_vl_ghia.append(float(data[0]))
			u_vl_ghia.append(float(data[row]))
	infile.close()
	
	infile = open(os.getcwd()+'/resources/ghia_et_al_1982/v_horizontal_line.dat', 'r')
	x_hl_ghia, v_hl_ghia = [], []
	index, row = 0, 0
	for line in infile:
		data = line.split()
		if index==0:
			for i in xrange(1, len(data)):
				if float(data[i]) == Solver.Re:
					row=i
				index += 1
		else:
			x_hl_ghia.append(float(data[0]))
			v_hl_ghia.append(float(data[row]))
	infile.close()


	y_vl = np.empty(Mesh.Ny, dtype=float)
	u_vl = np.empty(Mesh.Ny, dtype=float)
	x_hl = np.empty(Mesh.Nx, dtype=float)
	v_hl = np.empty(Mesh.Nx, dtype=float)
	I, J = 0, 0
	for i in xrange(Mesh.Nx*Mesh.Ny):
		if Mesh.x[i%Mesh.Nx] == Mesh.x[Mesh.Nx/2]:
			y_vl[I] = Mesh.y[i/Mesh.Nx]
			u_vl[I] = u.field[i]
			I += 1
		if Mesh.y[i/Mesh.Nx] == Mesh.y[Mesh.Ny/2]:
			x_hl[J] = Mesh.x[i%Mesh.Nx]
			v_hl[J] = v.field[i]
			J += 1

	plt.figure(num=None)
	plt.grid(True)
	plt.xlabel(r'$y$', fontsize=20)
	plt.ylabel(r'$u$', fontsize=20)
	plt.plot(y_vl, u_vl, 'b-', linewidth=2)
	plt.plot(y_vl_ghia, u_vl_ghia, 'ro', markersize=6)
	plt.legend(['pyIBM', 'Ghia et al. (1982)'], loc='best', prop={'size':16})
	plt.savefig(Case.path+'/images/velocity_vertical_line_'+str(Solver.start+Solver.nt)+'.png')
	plt.clf()
	plt.close()

	plt.figure(num=None)
	plt.grid(True)
	plt.xlabel(r'$x$', fontsize=20)
	plt.ylabel(r'$v$', fontsize=20)
	plt.plot(x_hl, v_hl, 'b-', linewidth=2)
	plt.plot(x_hl_ghia, v_hl_ghia, 'ro', markersize=6)
	plt.legend(['pyIBM', 'Ghia et al. (1982)'], loc='best', prop={'size':16})
	plt.savefig(Case.path+'/images/velocity_horizontal_line_'+str(Solver.start+Solver.nt)+'.png')
	plt.clf()
	plt.close()

if __name__ == '__main__':
	print '\n\t----- pyIBM - Cavity flow - Post-processing -----\n'
	main(sys.argv)
	print '\n\t----- pyIBM - END -----\n'
