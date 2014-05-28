# script: $pyIBM/scripts/post_processing_cavity.py
# Olivier Mesnard (mesnardo@gwu.edu)
# BarbaGroup (lorenabarba.com)

import os
import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0,'./src')
from case import *
from mesh import *
from solver import *
from variable import *


def main():
	"""Script to plot velocity at center lines
	and to compare with Ghia et al. (1982) experimental data.
	"""
	parser = argparse.ArgumentParser(description='Plots velocity at the centerline in both horizontal and vertical directions of the cavity')
	parser.add_argument('-p', '--path', dest='path', 
						help='path of the case folder', type=str)
	args = parser.parse_args()

	# creates the case
	Case(args.path)

	# generates the mesh
	mesh = Mesh()

	# initializes the solver
	Solver()

	# creates the velocity variables
	u = Variable('u')
	v = Variable('v')
	p = Variable('p')

	# reads variables' field at last iteration saved
	print '{Post-processing}: Comparison with Ghia et al. (1982)'
	Solver.ite = Solver.start + Solver.nt
	u.read()
	v.read()
	p.read()
	
	# computes velocity on the centered vertical line of the cavity
	file_path = os.getcwd()+'/resources/ghia_et_al_1982/u_vertical_line.dat'
	with open(file_path, 'r') as file_name:
		data = np.genfromtxt(file_name, 
							 dtype=float, names=True, unpack=True)
		y_vl_ghia = data['y']
		u_vl_ghia = data[str(Solver.Re)]
	
	# computes velocity on the centered horizontal line of the cavity
	file_path = os.getcwd()+'/resources/ghia_et_al_1982/v_horizontal_line.dat'
	with open(file_path, 'r') as file_name:
		data = np.genfromtxt(file_name, 
							 dtype=float, names=True, unpack=True)
		x_hl_ghia = data['x']
		v_hl_ghia = data[str(Solver.Re)]
	
	# creates the spatial arrays
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

	# creates the figure of the vertical line
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

	# creates the figure of the horizontal line
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
	main()
	print '\n\t----- pyIBM - END -----\n'
