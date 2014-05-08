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
			pmin = p.field.min()
			pmax = p.field.max()
			levels = np.arange(pmin,pmax+(pmax-pmin)/nc,(pmax-pmin)/nc)
			if (len(levels)>0):
				cont = plt.contourf(Mesh.x,Mesh.y,p.field.reshape(Mesh.Ny,Mesh.Nx),\
							levels,extend='both',cmap=cm.jet)
				cbar = plt.colorbar(cont)
				cbar.set_label(p.name)
			plt.streamplot(Mesh.x,Mesh.y,\
						u.field.reshape(Mesh.Ny,Mesh.Nx),\
						v.field.reshape(Mesh.Ny,Mesh.Nx))
			plt.xlim(Mesh.xmin,Mesh.xmax)
			plt.ylim(Mesh.ymin,Mesh.ymax)
			plt.title('pressure - '+str(ite))
			if (os.path.isdir('./images')==False):
				os.system('mkdir ./images')
			plt.savefig('./images/'+p.name+str('%05d'%(ite,))+'.png')
			plt.clf()
			plt.close()
	
	print '{Post-processing}: Comparison with Ghia et al. (1982)'
	u.read(Solver.start+Solver.nt)
	v.read(Solver.start+Solver.nt)
	p.read(Solver.start+Solver.nt)
	
	inFile = open('./ghia_et_al/u_verticalLine.dat','r')
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
	inFile = open('./ghia_et_al/v_horizontalLine.dat','r')
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

	plt.figure()
	plt.grid(True)
	plt.xlabel('y',fontsize=16)
	plt.ylabel('u',fontsize=16)
	plt.plot(y_VL,u_VL,'b-',linewidth=2)
	plt.plot(y_VL_ghia,u_VL_ghia,'ro',markersize=6)
	plt.legend(['present','Ghia et al. (1982)'],'best',prop={'size':16})
	plt.savefig('./images/velocityVL_'+str(Solver.start+Solver.nt)+'.png')
	plt.clf()
	plt.close()
	plt.figure()
	plt.grid(True)
	plt.xlabel('x',fontsize=16)
	plt.ylabel('v',fontsize=16)
	plt.plot(x_HL,v_HL,'b-',linewidth=2)
	plt.plot(x_HL_ghia,v_HL_ghia,'ro',markersize=6)
	plt.legend(['present','Ghia et al. (1982)'],'best',prop={'size':16})
	plt.savefig('./images/velocityHL_'+str(Solver.start+Solver.nt)+'.png')
	plt.clf()
	plt.close()

if (__name__=='__main__'):
	print '\n\t----- pyIBM - Cavity flow - Post-processing -----\n'
	main(sys.argv)
	print '\n\t----- pyIBM - END -----\n'
