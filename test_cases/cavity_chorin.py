# pyIBM - cavity flow - projection method of Chorin
# Olivier Mesnard
# mesnardo@gwu.edu

import numpy as np
import sys
import os
import matplotlib.pyplot as plt
from matplotlib import cm
import pyamg

sys.path.insert(0,'../src')
from mesh import *
from solver import *
from variable import *
from matrix import *
from operations import *
from poisson import *

import timeInfo as timeInfo

def main(arg):
	if (len(arg)==1): case = 'cavity_chorin'
	else: case = str(arg[1])


	print '\n{Meshing}'
	tic = timeInfo.start()
	mesh = Mesh('./'+case+'/_infoMesh.yaml')
	mesh.generate()
	timeInfo.stop(tic,'Meshing')
	print '{Writing mesh}'
	mesh.write()
	print '{Plotting mesh}'
	mesh.plot(isShow=False)

	print '\n{Creating solver}'
	solver = Solver('./'+case+'/_infoSolver.yaml')

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
	
	p = Variable('p')
	p.assembleMatrix('laplacian',scheme='central')
	p.assembleMatrix('gradientX',scheme='central',direction='x')
	p.assembleMatrix('gradientY',scheme='central',direction='y')
	timeInfo.stop(tic,'Assemble matrices')

	poissonP = Poisson(p,'./'+case+'/_infoSolver.yaml')

	u.write(Solver.start)
	v.write(Solver.start)
	p.write(Solver.start)

	b = np.empty(Mesh.Nx*Mesh.Ny,dtype=np.float64)

	if (Solver.scheme == 'Euler'): order=1
	elif (Solver.scheme == 'RK3'): order=3
	else: order=1
	alpha = [1.,3./4,1./3]
	beta = [0.,1./4,2./3]
	gamma = [1.,1./4,2./3]

	tic = timeInfo.start()
	while(Solver.ite<Solver.start+Solver.nt):
		Solver.ite += 1
		print '\nIteration ',Solver.ite,' - Time = ',Solver.ite*Solver.dt
		u.prev,v.prev,p.prev = u.field,v.field,p.field

		u.field[:] += Solver.dt*(
					1./Solver.Re*lap(u)[:]\
					-u.field[:]*grad(u,'x')[:]\
					-v.field[:]*grad(u,'y')[:])
		
		v.field[:] += Solver.dt*(
					1./Solver.Re*lap(v)[:]\
					-u.field[:]*grad(v,'x')[:]\
					-v.field[:]*grad(v,'y')[:])

		b[:] = 1./Solver.dt*(grad(u,'x')[:]+grad(v,'y')[:])

		p.field = poissonP.solve(p.laplacian.mat,b-p.laplacian.bcVect,p.field)
		print '{Poisson} Number of iterations: P=',poissonP.ite
			
		u.field[:] += -Solver.dt*grad(p,'x')[:]
		v.field[:] += -Solver.dt*grad(p,'y')[:]

		if (Solver.ite%Solver.writeEvery==0):
			print '\n{Writing results}'
			u.write(Solver.ite)
			v.write(Solver.ite)
			p.write(Solver.ite)

	timeInfo.stop(tic,'DONE')

if (__name__=='__main__'):
	print '\n\t----- pyIBM - Cavity flow - Step11 - CFD class -----\n'
	main(sys.argv)
	print '\n\t----- pyIBM - END -----\n'
