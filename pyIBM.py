#!/usr/env/bin python

##############################################
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

	# get case path
	pwd = os.getcwd()
	case_name = arg[1]
	case_path = pwd+'/'+case_name

	# create mesh
	print '\n{Meshing}'
	tic = timeInfo.start()
	mesh = Mesh(case_path+'/_infoMesh.yaml')
	mesh.generate()
	timeInfo.stop(tic, 'Meshing')
	print '{Writing mesh}'
	mesh.write()

	# create immersed boundary
	if Mesh.is_body:
		print '{Creating body}'
		body = Body(case_path+'/_infoBody.yaml')

	print '{Plotting mesh}'
	mesh.plot(body if Mesh.is_body else None, is_show=False)

	# create solver
	print '\n{Creating solver}'
	solver = Solver(case_path+'/_infoSolver.yaml')

	# create variables and related matrices
	print '\n{Assembling matrices}'
	tic = timeInfo.start()
	print 'u'
	u = Variable('u')
	u.assemble_matrix('laplacian', scheme='central')
	u.assemble_matrix('gradient_x', scheme='central', direction='x')
	u.assemble_matrix('gradient_y', scheme='central', direction='y')

	print 'v'
	v = Variable('v')
	v.assemble_matrix('laplacian', scheme='central')
	v.assemble_matrix('gradient_x', scheme='central', direction='x')
	v.assemble_matrix('gradient_y', scheme='central', direction='y')
	
	print 'p'
	p = Variable('p')
	p.assemble_matrix('laplacian', scheme='central')
	p.assemble_matrix('gradient_x', scheme='central', direction='x')
	p.assemble_matrix('gradient_y', scheme='central', direction='y')
	timeInfo.stop(tic, 'Assembling matrices')

	# create Poisson solver
	poisson_p = Poisson(p, case_path+'/_infoSolver.yaml')

	# initalization of the RHS of the Poisson equation
	b = np.empty(Mesh.Nx*Mesh.Ny, dtype=float)

	# open file to store force coefficients
	outfile = open(case_path+'/forceCoeffs.dat', ('w' if Solver.start == 0 else 'a'))

	tic = timeInfo.start()
	while Solver.ite < Solver.start + Solver.nt:
		Solver.ite += 1
		print '\nIteration ', Solver.ite, ' - Time = ', Solver.ite*Solver.dt
		
		u.field[:] += Solver.dt*(\
					+1./Solver.Re*lap(u)[:]\
					-u.field[:]*grad(u, 'x')[:]\
					-v.field[:]*grad(u, 'y')[:])
		v.field[:] += Solver.dt*(\
					+1./Solver.Re*lap(v)[:]\
					-u.field[:]*grad(v, 'x')[:]\
					-v.field[:]*grad(v, 'y')[:])
		
		# immersed boundary method
		if Mesh.is_body:
			ibm(body, u, v)
			outfile.write(str(Solver.ite*Solver.dt)+'\t'\
						+str(body.cl)+'\t'+str(body.cd)+'\n')
			print '{Body} \t Cl = %.3f \t Cd = %.3f' % (body.cl,body.cd)
		
		# solve the Poisson equation for pressure
		b[:] = 1./Solver.dt*(grad(u, 'x')[:]+grad(v, 'y')[:])
		p.field = poisson_p.solve(p.laplacian.mat, b-p.laplacian.bc_vect, p.field)
		
		print '{Poisson} Number of iterations: ', poisson_p.ite
		
		# update velocity field
		u.field[:] -= Solver.dt*grad(p, 'x')[:]
		v.field[:] -= Solver.dt*grad(p, 'y')[:]
		
		# write variable fields
		if Solver.ite%Solver.write_every == 0:
			print '\n{Writing results}'
			u.write()
			v.write()
			p.write()

	timeInfo.stop(tic, 'DONE')
	outfile.close()

if __name__ == '__main__':
	print '\n\t----- pyIBM - START -----\n'
	main(sys.argv)
	print '\n\t----- pyIBM - END -----\n'
