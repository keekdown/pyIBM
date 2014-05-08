# pyIBM - cavity flow - postProcessingCavity.py
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

import timeInfo as timeInfo

def getInfo(case):
	mesh = Mesh(case+'/_infoMesh.yaml')
	mesh.read()
	solver = Solver(case+'/_infoSolver.yaml')
	u = Variable('u')
	v = Variable('v')
	u.read(Solver.start+Solver.nt)
	v.read(Solver.start+Solver.nt)
	
	y_VL = np.empty(Mesh.Ny,dtype=np.float64)
	u_VL = np.empty(Mesh.Ny,dtype=np.float64)
	x_HL = np.empty(Mesh.Nx,dtype=np.float64)
	v_HL = np.empty(Mesh.Nx,dtype=np.float64)
	I,J = 0,0
	for i in range(Mesh.Nx*Mesh.Ny):
		if (Mesh.x[i%Mesh.Nx] == Mesh.x[Mesh.Nx/2]):
			y_VL[I] = Mesh.y[i/Mesh.Nx]
			u_VL[I] = u.field[i]
			I += 1
		if (Mesh.y[i/Mesh.Nx] == Mesh.y[Mesh.Ny/2]):
			x_HL[J] = Mesh.x[i%Mesh.Nx]
			v_HL[J] = v.field[i]
			J += 1
	return [y_VL,u_VL,x_HL,v_HL]

def main(arg):

	[y_VL_euler,u_VL_euler,x_HL_euler,v_HL_euler] = getInfo('cavity_Re400_Euler')
	[y_VL_rk3,u_VL_rk3,x_HL_rk3,v_HL_rk3] = getInfo('cavity_Re400_RK3')

	inFile = open('./step11_CFDclass/ghia_et_al/u_verticalLine.dat','r')
	y_VL_ghia,u_VL_ghia = [],[]
	index,row = 0,0
	for line in inFile:
		data = line.split()
		if (index==0):
			for i in range(1,len(data)):
				if (float(data[i])==Solver.Re): row=i
				index += 1
		else:
			y_VL_ghia.append(float(data[0]))
			u_VL_ghia.append(float(data[row]))
	inFile.close()
	inFile = open('./step11_CFDclass/ghia_et_al/v_horizontalLine.dat','r')
	x_HL_ghia,v_HL_ghia = [],[]
	index,row = 0,0
	for line in inFile:
		data = line.split()
		if (index==0):
			for i in range(1,len(data)):
				if (float(data[i])==Solver.Re): row=i
				index += 1
		else:
			x_HL_ghia.append(float(data[0]))
			v_HL_ghia.append(float(data[row]))
	inFile.close()

	plt.figure()
	plt.grid(True)
	plt.xlabel('y',fontsize=16)
	plt.ylabel('u',fontsize=16)
	plt.plot(y_VL_euler,u_VL_euler,'b-',linewidth=2)
	plt.plot(y_VL_rk3,u_VL_rk3,'g-',linewidth=2)
	plt.plot(y_VL_ghia,u_VL_ghia,'ro',markersize=6)
	plt.legend(['Euler','RK3','Ghia et al. (1982)'],'best',prop={'size':16})
	plt.savefig('./velocityVL_Re'+str(Solver.Re)+'_'+str(Solver.start+Solver.nt)+'.png')
	plt.clf()
	plt.close()
	plt.figure()
	plt.grid(True)
	plt.xlabel('x',fontsize=16)
	plt.ylabel('v',fontsize=16)
	plt.plot(x_HL_euler,v_HL_euler,'b-',linewidth=2)
	plt.plot(x_HL_rk3,v_HL_rk3,'g-',linewidth=2)
	plt.plot(x_HL_ghia,v_HL_ghia,'ro',markersize=6)
	plt.legend(['Euler','RK3','Ghia et al. (1982)'],'best',prop={'size':16})
	plt.savefig('./velocityHL_Re'+str(Solver.Re)+'_'+str(Solver.start+Solver.nt)+'.png')
	plt.clf()
	plt.close()

if (__name__=='__main__'):
	print '\n\t----- pyIBM - Cavity flow - Post-processing -----\n'
	main(sys.argv)
	print '\n\t----- pyIBM - END -----\n'
