# script: $pyIBM/scripts/post_processing.py
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
from body import *
from solver import *
from variable import *
from operations import *

def plot_pressure(p, body=None, limits=None):
	'''Plot pressure field on the mesh.'''

	# reads variable
	p.read()
	
	# initializes figure
	plt.figure(num=None)
	plt.grid(False)
	plt.xlabel(r'$x$', fontsize=20)
	plt.xlabel(r'$y$', fontsize=20)

	# creates contour of pressure
	nc = 100
	pmin, pmax = p.field.min(), p.field.max()
	pmin, pmax = -0.5, 0.5
	if pmin == pmax:
		pmin, pmax = -0.5, 0.5
	levels = np.linspace(pmin, pmax, nc)
	
	cont = plt.contourf(Mesh.x, Mesh.y, p.field.reshape(Mesh.Ny, Mesh.Nx),
						levels, extend='both', cmap=cm.jet)
	cbar = plt.colorbar(cont)
	cbar.set_label('pressure')

	# plots body if there is one
	if Mesh.is_body:
		plt.plot(np.append(body.x,body.x[0]), np.append(body.y, body.y[0]), 'k', ls='-', lw=1)

	# sets axis limits
	plt.xlim(limits[0], limits[1])
	plt.ylim(limits[2], limits[3])

	# inserts title and saves figure
	plt.title('pressure - '+str(Solver.ite))
	plt.savefig(Case.path+'/images/'+'pressure'+str('%04d'%(Solver.ite,))+'.png')
	
	plt.clf()
	plt.close()

def plot_velocity(u, v, body=None, limits=None):
	'''Plot velocity field on the mesh.'''
	
	# reads variables
	u.read()
	v.read()
	
	# initializes figure
	plt.figure(num=None)
	plt.grid(False)
	plt.xlabel(r'$x$', fontsize=20)
	plt.xlabel(r'$y$', fontsize=20)

	# creates contour of velocity magnitude and adds streams
	nc = 100
	magn = np.sqrt(u.field**2+v.field**2)
	umin, umax = magn.min(), magn.max()
	umin, umax = 0., 1.
	levels = np.linspace(umin, umax, nc)

	cont = plt.contourf(Mesh.x, Mesh.y, magn.reshape(Mesh.Nx, Mesh.Ny),
						levels, extend='both', cmap=cm.jet)
	cbar = plt.colorbar(cont)
	cbar.set_label('velocity')

	plt.streamplot(Mesh.x, Mesh.y,
				   u.field.reshape(Mesh.Ny, Mesh.Nx),
				   v.field.reshape(Mesh.Ny, Mesh.Nx))

	# plots body if there is one
	if Mesh.is_body:
		plt.plot(np.append(body.x,body.x[0]), np.append(body.y, body.y[0]), 'k', ls='-', lw=1)

	# sets axis limits
	plt.xlim(limits[0], limits[1])
	plt.ylim(limits[2], limits[3])

	# inserts title and saves the figure
	plt.title('velocity - '+str(Solver.ite))
	plt.savefig(Case.path+'/images/'+'velocity'+str('%04d'%(Solver.ite,))+'.png')
	
	plt.clf()
	plt.close()


def plot_vorticity(u, v, body=None, limits=None):
	'''Plot vorticity field on the mesh.'''
	
	# reads the variables
	u.read()
	v.read()

	# computes vorticity
	w = grad(v,'x') - grad(u,'y')

	# initializes figure
	plt.figure(num=None)
	plt.grid(False)
	plt.xlabel(r'$x$', fontsize=20)
	plt.xlabel(r'$y$', fontsize=20)

	# creates contour of vorticity
	nc = 100
	wmin, wmax = w.min(), w.max()
	wmin, wmax = -1.0, 1.0
	if wmin == wmax:
		wmin, wmax = -1.0, 1.0
	levels = np.linspace(wmin, wmax, nc)
	
	cont = plt.contourf(Mesh.x, Mesh.y, w.reshape(Mesh.Ny, Mesh.Nx),
						levels, extend='both', cmap=cm.jet)
	cbar = plt.colorbar(cont)
	cbar.set_label('vorticity')

	# plots body of there is one
	if Mesh.is_body:
		plt.plot(np.append(body.x,body.x[0]), np.append(body.y, body.y[0]), 'k', ls='-', lw=1)

	# sets axis limits
	plt.xlim(limits[0], limits[1])
	plt.ylim(limits[2], limits[3])

	# inserts title and saves figure
	plt.title('vorticity - '+str(Solver.ite))
	plt.savefig(Case.path+'/images/'+'vorticity'+str('%04d'%(Solver.ite,))+'.png')
	
	plt.clf()
	plt.close()


def main():
	'''Plot either pressure, velocity or vorticity,
	at every time saved in the case folder.
	'''

	# list of command-line arguments
	parser = argparse.ArgumentParser(description='Plots pressure, velocity and/or vorticity')
	parser.add_argument('-p', '--path', dest='path', 
						help='path of the case folder', type=str)
	parser.add_argument('-v', '--variable', dest='variable', 
						help='list of variables to plot (pressure, velocity, vorticity)', 
						nargs='+', type=str, default=['pressure', 'velocity', 'vorticity'])
	parser.add_argument('-z', '--zoom', dest='zoom', 
						help='limits of the plot', nargs='+', type=float)
	parser.add_argument('-t', '--time', dest='time', 
						help='time-step to plot / List: (min, max, every) or specific time', 
						nargs='+', type=int)
	args = parser.parse_args()

	# creates the case
	Case(args.path)
	
	# generates the mesh
	mesh = Mesh()
	
	# zoom on a given window
	if args.zoom != None:
		limits = args.zoom
	else:
		limits = [Mesh.xmin, Mesh.xmax, Mesh.ymin, Mesh.ymax]

	# creates potential bodies
	body = Body() if Mesh.is_body else None

	# initializes the solver
	Solver()
	
	# changes timesteps to plot if argument specified
	if args.time != None:
		if len(args.time) == 1:
			Solver.start, Solver.nt = args.time[0], 0
		else:
			Solver.start = args.time[0]
			Solver.write_every = args.time[2]
			Solver.nt = args.time[1] - Solver.start + Solver.write_every
			Solver.ite = Solver.start - Solver.write_every

	# chooses variables to plot if argument specified
	if args.variable != None:
		variables = args.variable
	else:
		variables = ['pressure', 'velocity', 'vorticity']

	print variables	

	# creates variables and necessary matrices for vorticity
	if 'pressure' in variables:
		p = Variable('p', skip_assemble=True)
	if 'velocity' in variables or 'vorticity' in variables:
		u = Variable('u', skip_assemble=True)
		v = Variable('v', skip_assemble=True)
		if 'vorticity' in variables:
			u.assemble_matrix('gradient', scheme='central', direction='y')
			v.assemble_matrix('gradient', scheme='central', direction='x')

	# create an 'images' folder in case folder if necessary
	if not os.path.isdir(Case.path+'/images'):
		os.system('mkdir '+Case.path+'/images')

	# loops over time
	for ite in xrange(Solver.start, Solver.start+Solver.nt, Solver.write_every):
		
		Solver.ite += Solver.write_every

		print 'Iteration ', Solver.ite

		# plots different variables
		if 'pressure' in variables:
			plot_pressure(p, body, limits)
		if 'velocity' in variables:
			plot_velocity(u, v, body, limits)
		if 'vorticity' in variables:
			plot_vorticity(u, v, body, limits)
		

if __name__ == '__main__':
	print '\n\t----- pyIBM - POST-PROCESSING -----\n'
	main()
	print '\n\t----- pyIBM - END -----\n'
