# pyIBM - cavity flow - postProcessingCavity.py
# Olivier Mesnard
# mesnardo@gwu.edu

import numpy as np
import sys
import os
import matplotlib.pyplot as plt
from matplotlib import cm

sys.path.insert(0,'../../src')
from mesh import *
from solver import *
from variable import *

import timeInfo as timeInfo

def main(arg):

	mesh = Mesh('./_infoMesh.yaml')
	mesh.read()
	mesh.plot()

	solver = Solver('./_infoSolver.yaml')

	u = Variable('u')
	v = Variable('v')
	p = Variable('p')

	start = Solver.start+Solver.writeEvery
	end = Solver.start+Solver.nt+Solver.writeEvery
	every = Solver.writeEvery
	for ite in range(start,end,every):
			print '{Post-processing}: Iteration ',ite

			u.read(ite)
			v.read(ite)
			p.read(ite)
			
			plt.figure()
			plt.grid(False)
			plt.xlabel('x',fontsize=16)
			plt.ylabel('y',fontsize=16)
			nc = 100
			magn = np.sqrt(u.field**2+v.field**2)
			umin = magn.min()
			umax = magn.max()
			levels = np.arange(umin,umax+(umax-umin)/nc,(umax-umin)/nc)
			if (len(levels)>0):
				cont = plt.contourf(Mesh.x,Mesh.y,magn.reshape(Mesh.Ny,Mesh.Nx),\
							levels,extend='both',cmap=cm.jet)
				cbar = plt.colorbar(cont)
				cbar.set_label('U')
			plt.streamplot(Mesh.x,Mesh.y,\
						u.field.reshape(Mesh.Ny,Mesh.Nx),\
						v.field.reshape(Mesh.Ny,Mesh.Nx))
			plt.xlim(Mesh.xmin,Mesh.xmax)
			plt.ylim(Mesh.ymin,Mesh.ymax)
			plt.title('velocity - '+str(ite))
			if (os.path.isdir('./images')==False):
				os.system('mkdir ./images')
			plt.savefig('./images/'+'velocity'+str('%05d'%(ite,))+'.png')
			plt.clf()
			plt.close()

if (__name__=='__main__'):
	print '\n\t----- pyIBM - Channel flow - Post-processing -----\n'
	main(sys.argv)
	print '\n\t----- pyIBM - END -----\n'
