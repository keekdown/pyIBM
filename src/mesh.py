# source: $pyIBM/src/mesh.py
# Olivier Mesnard (mesnardo@gwu.edu)
# BarbaGroup (lorenabarba.com)
		
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import os
import sys
import yaml		

class Subdomain:
	def __init__(self,info_subdomain):
		self.end = info_subdomain['end']
		self.N = info_subdomain['cells']
		self.gamma = info_subdomain['stretchRatio']

class Direction:
	def __init__(self,info_direction):
		self.name = info_direction['direction']
		self.start = info_direction['start']
		self.end = info_direction['subdomains'][-1]['end']
		self.subdomain = []
		for i in range(len(info_direction['subdomains'])):
			self.subdomain.append(Subdomain(info_direction['subdomains'][i]))
		self.N = sum(self.subdomain[i].N for i in range(len(self.subdomain)))

	def generate(self):
		# WARNING: only work for uniform mesh
		# Adaptation to stretched mesh will be done later
		self.coord = np.empty(self.N,dtype=float)
		self.delta = np.empty(self.N,dtype=float)
		h = (self.end-self.start)/(self.N-1)
		for i in range(self.N):
			self.delta[i] = h
			self.coord[i] = self.start + i*h


class Mesh:
	def __init__(self,info_mesh):
		Mesh.info_mesh = info_mesh
		Mesh.case_path = os.path.dirname(os.path.abspath(self.info_mesh))
		Mesh.is_body = False
		infile = open(self.info_mesh,'r')
		info = yaml.load(infile)
		infile.close()
		Mesh.direction = []
		for i in range(len(info)):
			if ('direction' in info[i]):
				self.direction.append(Direction(info[i]))
			elif ('body' in info[i]):
				Mesh.is_body = True
		Mesh.xmin,Mesh.xmax = Mesh.direction[0].start,Mesh.direction[0].end
		Mesh.ymin,Mesh.ymax = Mesh.direction[1].start,Mesh.direction[1].end
		Mesh.Nx = Mesh.direction[0].N
		Mesh.Ny = Mesh.direction[1].N
		print '-> Number of points in the x-direction: ',Mesh.Nx
		print '-> Number of points in the y-direction: ',Mesh.Ny

	def generate(self):
		for item in Mesh.direction:
			item.generate()
		Mesh.x,Mesh.y = Mesh.direction[0].coord,Mesh.direction[1].coord
		Mesh.dx,Mesh.dy = Mesh.direction[1].delta,Mesh.direction[1].delta

	def write(self):
		outfile = open(Mesh.case_path+'/mesh.dat','w')
		outfile.write(str(Mesh.Nx)+'\t'+str(Mesh.Ny)+'\n')
		for i in range(Mesh.Nx):
			outfile.write(str(Mesh.x[i])+'\t'+str(Mesh.dx[i])+'\n')
		for i in range(Mesh.Ny):
			outfile.write(str(Mesh.y[i])+'\t'+str(Mesh.dy[i])+'\n')
		outfile.close()

	def read(self):
		infile = open(Mesh.case_path+'/mesh.dat','r')
		index,i = 0,0
		for line in infile:
			data = line.split()
			if (index == 0):
				Mesh.Nx,Mesh.Ny = int(data[0]),int(data[1])
				Mesh.x = np.empty(Mesh.Nx,dtype=np.float64)
				Mesh.y = np.empty(Mesh.Ny,dtype=np.float64)
				Mesh.dx = np.empty(Mesh.Nx,dtype=np.float64)
				Mesh.dy = np.empty(Mesh.Ny,dtype=np.float64)
				index += 1
			else:
				if (i<Mesh.Nx):
					Mesh.x[i],Mesh.dx[i] = float(data[0]),float(data[1])
				else:
					Mesh.y[i-Mesh.Nx],Mesh.dy[i-Mesh.Nx] = float(data[0]),float(data[1])
				i += 1
		infile.close()

	def plot(self,body=None,is_show=False):
		plt.figure(num=None)
		plt.grid(False)
		plt.xlabel('x',fontsize=16)
		plt.ylabel('y',fontsize=16)
		for i in range(Mesh.Nx):
			plt.axvline(Mesh.x[i],Mesh.ymin,Mesh.ymax,linewidth=1)
		for i in range(Mesh.Ny):
			plt.axhline(Mesh.y[i],Mesh.xmin,Mesh.xmax,linewidth=1)
		plt.xlim(Mesh.xmin,Mesh.xmax)
		plt.ylim(Mesh.ymin,Mesh.ymax)
		if (body != None):
			plt.plot(body.xk,body.yk,'ko-',lw=1,markersize=4)
			for k in range(body.Nk):
				plt.plot(Mesh.x[body.neighbor[k]%Mesh.Nx],Mesh.y[body.neighbor[k]/Mesh.Nx],\
						'ro',markersize=4)
			plt.title('MESH: '+str(Mesh.Nx)+'x'+str(Mesh.Ny)+' / IB: '+str(body.Nk))
			plt.savefig(Mesh.case_path+'/mesh_'+str(Mesh.Nx)+'_'+str(Mesh.Ny)+'_'+str(body.Nk)+'.png')
		else:
			plt.title('MESH: '+str(Mesh.Nx)+'x'+str(Mesh.Ny))
			if not(os.path.isdir(Mesh.case_path+'/images')):
				os.system('mkdir '+Mesh.case_path+'/images')
			plt.savefig(Mesh.case_path+'/images/mesh_'+str(Mesh.Nx)+'_'+str(Mesh.Ny)+'.png')
		if (is_show):
			plt.show()
		plt.clf()
		plt.close()
