# pyIBM - step11_CFDclass.py - 2D inviscid Burgers' equations
# Olivier Mesnard
# mesnardo@gwu.edu

import numpy as np
import sys
import os
import matplotlib.pyplot as plt
from matplotlib import cm

sys.path.insert(0,'../src')
from mesh import *
from solver import *
from variable import *
from matrix import *
from operations import *
from poisson import *

import timeInfo as timeInfo

def main(arg):

	print '\n{Meshing}'
	tic = timeInfo.start()
	mesh = Mesh('./step8_CFDclass/_infoMesh.yaml')
	mesh.generate()
	timeInfo.stop(tic,'Meshing')
	print '{Writing mesh}'
	mesh.write()
	print '{Plotting mesh}'
	mesh.plot()

	print '\n{Creating solver}'
	solver = Solver('./step8_CFDclass/_infoSolver.yaml')

	print '\n{Assembling matrices}'
	tic = timeInfo.start()
	u = Variable('u')
	u.assembleMatrix('laplacian',scheme='central')
	u.assembleMatrix('gradientX',scheme='backward',direction='x')
	u.assembleMatrix('gradientY',scheme='backward',direction='y')

	v = Variable('v')
	v.assembleMatrix('laplacian',scheme='central')
	v.assembleMatrix('gradientX',scheme='backward',direction='x')
	v.assembleMatrix('gradientY',scheme='backward',direction='y')
	
	timeInfo.stop(tic,'Assemble matrices')

	# Initial Conditions
	for j in range(Mesh.Ny):
		for i in range(Mesh.Nx):
			if (0.5<=Mesh.x[i]<=1.0 and 0.5<=Mesh.y[j]<=1.0):
				u.field[j*Mesh.Nx+i] = 2.0
				v.field[j*Mesh.Nx+i] = 2.0

	u.write(Solver.start)
	v.write(Solver.start)

	tic = timeInfo.start()
	while(Solver.ite<Solver.start+Solver.nt):
		Solver.ite += 1
		print '\nIteration ',Solver.ite,' - Time = ',Solver.ite*Solver.dt
		u.prev,v.prev = u.field,v.field

		u.field[:] = u.prev[:] + Solver.dt*(\
					+ 1./Solver.Re*lap(u)[:]\
					- u.prev[:]*grad(u,'x')[:]\
					- v.prev[:]*grad(u,'y')[:])
		v.field[:] = v.prev[:] + Solver.dt*(\
					+ 1./Solver.Re*lap(v)[:]\
					- u.prev[:]*grad(v,'x')[:]\
					- v.prev[:]*grad(v,'y')[:])
			
		if (Solver.ite%Solver.writeEvery==0):
			print '\n{Writing results}'
			u.write(Solver.ite)
			v.write(Solver.ite)

			plt.figure()
			plt.grid(True)
			plt.xlabel('x',fontsize=16)
			plt.ylabel('y',fontsize=16)
			nc = 100
			U = np.sqrt(u.field**2+v.field**2)
			Umin = U.min()
			Umax = U.max()
			levels = np.arange(Umin,Umax+(Umax-Umin)/nc,(Umax-Umin)/nc)
			if (len(levels)>0):
				cont = plt.contourf(Mesh.x,Mesh.y,U.reshape(Mesh.Ny,Mesh.Nx),\
							levels,extend='both',cmap=cm.jet)
				cbar = plt.colorbar(cont)
				cbar.set_label('U')
			plt.xlim(Mesh.xmin,Mesh.xmax)
			plt.ylim(Mesh.ymin,Mesh.ymax)
			plt.title('velocity magnitude - '+str(Solver.ite))
			if (os.path.isdir('./step8_CFDclass/images')==False):
				os.system('mkdir ./step8_CFDclass/images')
			plt.savefig('./step8_CFDclass/images/U'+str('%04d'%(Solver.ite,))+'.png')
			plt.clf()
			plt.close()

	timeInfo.stop(tic,'DONE')

if (__name__=='__main__'):
	print '\n\t----- pyIBM - 2D Inviscid Burgers - Step8 - CFD class -----\n'
	main(sys.argv)
	print '\n\t----- pyIBM - END -----\n'
