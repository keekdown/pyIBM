# pyIBM  - test case: poisson solver - step10 from CFD class
# Olivier Mesnard
# mesnardo@gwu.edu

import numpy as np
import sys
import os
import matplotlib.pyplot as plt

sys.path.insert(0,'../src')
from mesh import *
from solver import *
from variable import *
from matrix import *
from poisson import *

import timeInfo as timeInfo

def plotResults(mesh,p):
	plt.figure()
	plt.grid(True)
	plt.xlabel('x',fontsize=24)
	plt.ylabel('y',fontsize=24)
	pmin,pmax = p.field.min(),p.field.max()
	nc = 100
	levels = np.arange(pmin,pmax+(pmax-pmin)/nc,(pmax-pmin)/nc)
	cont = plt.contourf(mesh.x,mesh.y,p.field.reshape(mesh.Ny,mesh.Nx),\
						levels,extend='both',cmap=cm.jet)
	cbar = plt.colorbar(cont)
	plt.xlim(mesh.xmin,mesh.xmax)
	plt.ylim(mesh.ymin,mesh.ymax)
	plt.title('pressure')
	plt.savefig('./step10_CFDclass/pressure.png')
	plt.clf()
	plt.close()

def main(arg):

	mesh = Mesh('./step10_CFDclass/_infoMesh.yaml')
	mesh.generate()
	mesh.write()
	mesh.plot()

	solver = Solver('./step10_CFDclass/_infoSolver.yaml')

	p = Variable('p')
	p.assembleMatrix('laplacian',scheme='central')

	poisson = Poisson('cg',1.0E-06,1000)

	b = np.zeros(Mesh.Nx*Mesh.Ny,dtype=np.float64)
	b[Mesh.Ny/4*Mesh.Nx+Mesh.Ny/4] = 100.
	b[Mesh.Ny/4*Mesh.Nx+3*Mesh.Ny/4] = 100.
	b[3*Mesh.Ny/4*Mesh.Nx+Mesh.Ny/4] = 100.
	b[3*Mesh.Ny/4*Mesh.Nx+3*Mesh.Ny/4] = 100.
	b[Mesh.Ny/4*Mesh.Nx+Mesh.Ny/2] = -100.
	b[Mesh.Ny/2*Mesh.Nx+Mesh.Ny/4] = -100.
	b[Mesh.Ny/2*Mesh.Nx+3*Mesh.Ny/4] = -100.
	b[3*Mesh.Ny/4*Mesh.Nx+Mesh.Ny/2] = -100.

	p.field = poisson.solve(p.laplacian.mat,b-p.laplacian.bcVect,p.field)
	print '\n{Poisson}: Number of iterations: ',poisson.ite

	plotResults(mesh,p)

	if ('--show' in arg): plt.show()

if (__name__=='__main__'):
	print '\n\t----- Poisson Equation -----\n'
	main(sys.argv)
	print '\n\t----- END -----\n'
