# script: $pyIBM/scripts/post_processing_cavity.py
# Olivier Mesnard (mesnardo@gwu.edu)
# BarbaGroup (lorenabarba.com)


import os
import sys
import argparse

import numpy as np
from matplotlib import pyplot as plt

sys.path.insert(0,'./src')
from case import Case
from mesh import Mesh
from parameters import Parameters
from solver import Solver
from input_output import read_inputs


def main():
	"""Script to plot velocity at center lines
	and to compare with Ghia et al. (1982) experimental data.
	"""
	# parse the command-line
	args = read_inputs()

	# create the case
	Case(args.path)

	# generate the mesh
	mesh = Mesh()

	# initialize the solver
	Solver(skip_assemble=True, skip_poisson=True)

	# read variables field at last iteration saved
	print '{Post-processing}: Comparison with Ghia et al. (1982)'
	Parameters.ite = Parameters.start + Parameters.nt
	Solver.u.read()
	Solver.v.read()
	Solver.p.read()
	
	# read velocity on the centered vertical line of the cavity
	file_path = os.getcwd()+'/resources/ghia_et_al_1982/u_vertical_line.dat'
	with open(file_path, 'r') as file_name:
		data = np.genfromtxt(file_name, 
							 dtype=float, names=True, unpack=True)
		y_vl_ghia = data['y']
		u_vl_ghia = data[str(Parameters.Re)]
	
	# read velocity on the centered horizontal line of the cavity
	file_path = os.getcwd()+'/resources/ghia_et_al_1982/v_horizontal_line.dat'
	with open(file_path, 'r') as file_name:
		data = np.genfromtxt(file_name, 
							 dtype=float, names=True, unpack=True)
		x_hl_ghia = data['x']
		v_hl_ghia = data[str(Parameters.Re)]
	
	# create arrays to store the solution
	y_vl = np.empty(Mesh.Ny, dtype=float)
	u_vl = np.empty(Mesh.Ny, dtype=float)
	x_hl = np.empty(Mesh.Nx, dtype=float)
	v_hl = np.empty(Mesh.Nx, dtype=float)
	I, J = 0, 0
	for i in xrange(Mesh.Nx*Mesh.Ny):
		if Mesh.x[i%Mesh.Nx] == Mesh.x[Mesh.Nx/2]:
			y_vl[I] = Mesh.y[i/Mesh.Nx]
			u_vl[I] = Solver.u.field[i]
			I += 1
		if Mesh.y[i/Mesh.Nx] == Mesh.y[Mesh.Ny/2]:
			x_hl[J] = Mesh.x[i%Mesh.Nx]
			v_hl[J] = Solver.v.field[i]
			J += 1

	# plot velcoity on the vertical line
	plt.figure(num=None)
	plt.grid(True)
	plt.xlabel(r'$y$', fontsize=18)
	plt.ylabel(r'$u$', fontsize=18)
	plt.plot(y_vl, u_vl, color='b', ls='-', lw=2.0)
	plt.scatter(y_vl_ghia, u_vl_ghia, color='r', marker='o', s=6)
	plt.legend(['pyIBM', 'Ghia et al. (1982)'], loc='best', prop={'size':16})
	plt.savefig(Case.images+'/velocity_vertical_line_%d.png' % (Parameters.start+Parameters.nt))
	plt.clf()
	plt.close()

	# plot velocity on the horizontal line
	plt.figure(num=None)
	plt.grid(True)
	plt.xlabel(r'$x$', fontsize=18)
	plt.ylabel(r'$v$', fontsize=18)
	plt.plot(x_hl, v_hl, color='b', ls='-', lw=2.0)
	plt.scatter(x_hl_ghia, v_hl_ghia, color='r', marker='o', s=6)
	plt.legend(['pyIBM', 'Ghia et al. (1982)'], loc='best', prop={'size':16})
	plt.savefig(Case.images+'/velocity_horizontal_line_%d.png' % (Parameters.start+Parameters.nt))
	plt.clf()
	plt.close()

if __name__ == '__main__':
	print '\n\t----- pyIBM - Cavity flow - Post-processing -----\n'
	main()
	print '\n\t----- pyIBM - END -----\n'
