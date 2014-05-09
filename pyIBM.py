#!/usr/env/bin python

##############################################
# source: $pyIBM/src/pyIBM.py                #
# pyIBM - Immersed Boundary Method in Python #
# Olivier Mesnard (mesnardo@gwu.edu)         #
# BarbaGroup (lorenabarba.com)               #
##############################################

import os
import sys
import numpy as np

sys.path.insert(0,'./src')
from mesh import *
from body import *
from solver import *
from variable import *
from matrix import *
from operations import *
from ibm import *
from poisson import *

import timeInfo as timeInfo

def main(arg):

	pwd = os.getcwd()
	case_name = arg[1]
	case_path = pwd+'/'+case_name

	print '\n{Meshing}'
	tic = timeInfo.start()
	mesh = Mesh(case_path+'/_infoMesh.yaml')
	mesh.generate()
	timeInfo.stop(tic,'Meshing')
	print '{Writing mesh}'
	mesh.write()

	body = None
	if (Mesh.is_body):
		print '{Creating body}'
		body = Body(case_path+'/_infoBody.yaml',Mesh)

	print '{Plotting mesh}'
	mesh.plot(body,is_show=False)

	print '\n{Creating solver}'
	solver = Solver(case_path+'/_infoSolver.yaml')

	print '\n{Assembling matrices}'
	tic = timeInfo.start()
	print 'u'
	u = Variable('u')
	u.assemble_matrix('laplacian',scheme='central')
	u.assemble_matrix('gradient_x',scheme='central',direction='x')
	u.assemble_matrix('gradient_y',scheme='central',direction='y')

	print 'v'
	v = Variable('v')
	v.assemble_matrix('laplacian',scheme='central')
	v.assemble_matrix('gradient_x',scheme='central',direction='x')
	v.assemble_matrix('gradient_y',scheme='central',direction='y')
	
	print 'p'
	p = Variable('p')
	p.assemble_matrix('laplacian',scheme='central')
	p.assemble_matrix('gradient_x',scheme='central',direction='x')
	p.assemble_matrix('gradient_y',scheme='central',direction='y')
	timeInfo.stop(tic,'Assembling matrices')

	poisson_p = Poisson(p,case_path+'/_infoSolver.yaml')

	u.write()
	v.write()
	p.write()

	b = np.empty(Mesh.Nx*Mesh.Ny,dtype=float)

	if (Solver.scheme == 'Euler'):
		order=1
	elif (Solver.scheme == 'RK3'):
		order=3
	else:
		order=1
	alpha = [1.,3./4,1./3]
	beta = [0.,1./4,2./3]
	gamma = [1.,1./4,2./3]

	if (Solver.start == 0):
		outfile = open(case_path+'/forceCoeffs.dat','w')
	else:
		outfile = open(case_path+'/forceCoeffs.dat','a')

	tic = timeInfo.start()
	while(Solver.ite<Solver.start+Solver.nt):
		Solver.ite += 1
		print '\nIteration ',Solver.ite,' - Time = ',Solver.ite*Solver.dt
		u.prev,v.prev,p.prev = u.field,v.field,p.field

		u.field[:] += Solver.dt*(\
					+1./Solver.Re*lap(u)[:]\
					-u.field[:]*grad(u,'x')[:]\
					-v.field[:]*grad(u,'y')[:])
		v.field[:] += Solver.dt*(\
					+1./Solver.Re*lap(v)[:]\
					-u.field[:]*grad(v,'x')[:]\
					-v.field[:]*grad(v,'y')[:])
		
		if (Mesh.is_body):
			fx = np.zeros(Mesh.Nx*Mesh.Ny,dtype=float)
			fy = np.zeros(Mesh.Nx*Mesh.Ny,dtype=float)
			if (body.is_moving):
				body.kinematics(Mesh)
			Cl,Cd = 0,0
			for i in range(1):
				body.u = interpolation(u.field,body)
				body.v = interpolation(v.field,body)
				body.fx[:] = (body.ud[:]-body.u[:])/Solver.dt
				body.fy[:] = (body.vd[:]-body.v[:])/Solver.dt
				for k in range(body.N):
					Cd += -2.*body.fx[k]\
							*Mesh.dx[body.neighbor[k]%Mesh.Nx]\
							*Mesh.dy[body.neighbor[k]/Mesh.Ny]
					Cl += -2.*body.fy[k]\
							*Mesh.dx[body.neighbor[k]%Mesh.Nx]\
							*Mesh.dy[body.neighbor[k]/Mesh.Nx]
				fx += distribution(body.fx,body)
				fy += distribution(body.fy,body)
			u.field[:] += Solver.dt*fx[:]
			v.field[:] += Solver.dt*fy[:]
			print 'Cl = ',Cl,'\tCd = ',Cd
			outfile.write(str(Solver.ite*Solver.dt)+'\t'\
						+str(Cl)+'\t'+str(Cd)+'\n')

		b[:] = 1./Solver.dt*(grad(u,'x')[:]+grad(v,'y')[:])
		p.field = poisson_p.solve(p.laplacian.mat,b-p.laplacian.bc_vect,p.field)
		print '{Poisson} Number of iterations: ',poisson_p.ite
		
		u.field[:] += -Solver.dt*grad(p,'x')[:]
		v.field[:] += -Solver.dt*grad(p,'y')[:]
		
		if (Solver.ite%Solver.write_every==0):
			print '\n{Writing results}'
			u.write()
			v.write()
			p.write()

	timeInfo.stop(tic,'DONE')
	outfile.close()

if (__name__=='__main__'):
	print '\n\t----- pyIBM - START -----\n'
	main(sys.argv)
	print '\n\t----- pyIBM - END -----\n'
