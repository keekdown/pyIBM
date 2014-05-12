# source: $pyIBM/src/mesh.py
# Olivier Mesnard (mesnardo@gwu.edu)
# BarbaGroup (lorenabarba.com)

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import yaml		

from case import Case

class Subdomain:
	def __init__(self, info_subdomain):
		self.end = info_subdomain['end']
		self.N = info_subdomain['cells']
		self.gamma = info_subdomain['stretchRatio']

class Direction:
	def __init__(self, info_direction):
		self.name = info_direction['direction']
		self.start = info_direction['start']
		self.end = info_direction['subdomains'][-1]['end']
		self.subdomain = []
		for i in xrange(len(info_direction['subdomains'])):
			self.subdomain.append(Subdomain(info_direction['subdomains'][i]))
		self.N = sum(self.subdomain[i].N for i in xrange(len(self.subdomain)))

	def generate(self):
		# WARNING: only work for uniform mesh
		# Adaptation to stretched mesh will be done later
		self.coord = np.empty(self.N, dtype=float)
		self.delta = np.empty(self.N, dtype=float)
		h = (self.end-self.start)/(self.N-1)
		for i in xrange(self.N):
			self.delta[i] = h
			self.coord[i] = self.start + i*h


class Mesh:
	def __init__(self):
		Mesh.is_body = False
		infile = open(Case.path+'/_infoMesh.yaml', 'r')
		info = yaml.load(infile)
		infile.close()
		Mesh.direction = []
		for i in xrange(len(info)):
			if 'direction' in info[i]:
				Mesh.direction.append(Direction(info[i]))
			elif 'body' in info[i]:
				Mesh.is_body = True
		Mesh.xmin, Mesh.xmax = Mesh.direction[0].start, Mesh.direction[0].end
		Mesh.ymin, Mesh.ymax = Mesh.direction[1].start, Mesh.direction[1].end
		Mesh.Nx = Mesh.direction[0].N
		Mesh.Ny = Mesh.direction[1].N

		self.generate()
		self.write()

		print '\n'
		print '-> Number of points in the x-direction: ', Mesh.Nx
		print '-> Number of points in the y-direction: ', Mesh.Ny
		print '\n'

	def generate(self):
		for d in Mesh.direction:
			d.generate()
		Mesh.x, Mesh.y = Mesh.direction[0].coord,Mesh.direction[1].coord
		Mesh.dx, Mesh.dy = Mesh.direction[1].delta,Mesh.direction[1].delta

	def write(self):
		outfile = open(Case.path+'/mesh.dat', 'w')
		outfile.write(str(Mesh.Nx)+'\t'+str(Mesh.Ny)+'\n')
		for i in xrange(Mesh.Nx):
			outfile.write(str(Mesh.x[i])+'\t'+str(Mesh.dx[i])+'\n')
		for i in xrange(Mesh.Ny):
			outfile.write(str(Mesh.y[i])+'\t'+str(Mesh.dy[i])+'\n')
		outfile.close()

	def read(self):
		infile = open(Case.path+'/mesh.dat', 'r')
		index, i = 0, 0
		for line in infile:
			data = line.split()
			if index == 0:
				Mesh.Nx, Mesh.Ny = int(data[0]), int(data[1])
				Mesh.x = np.empty(Mesh.Nx, dtype=float)
				Mesh.y = np.empty(Mesh.Ny, dtype=float)
				Mesh.dx = np.empty(Mesh.Nx, dtype=float)
				Mesh.dy = np.empty(Mesh.Ny, dtype=float)
				index += 1
			else:
				if i < Mesh.Nx:
					Mesh.x[i], Mesh.dx[i] = float(data[0]), float(data[1])
				else:
					Mesh.y[i-Mesh.Nx], Mesh.dy[i-Mesh.Nx] = float(data[0]), float(data[1])
				i += 1
		infile.close()

	def plot(self,body=None, is_show=False):
		if not os.path.isdir(Case.path+'/images'):
			os.system('mkdir '+Case.path+'/images')
		plt.figure(num=None)
		plt.grid(False)
		plt.xlabel('x', fontsize=16)
		plt.ylabel('y', fontsize=16)
		for i in xrange(Mesh.Nx):
			plt.axvline(Mesh.x[i])
		for i in xrange(Mesh.Ny):
			plt.axhline(Mesh.y[i])
		plt.xlim(Mesh.xmin, Mesh.xmax)
		plt.ylim(Mesh.ymin, Mesh.ymax)
		if body != None:
			plt.plot(np.append(body.x, body.x[0]), np.append(body.y, body.y[0]), 'ko-', lw=1, markersize=4)
			for k in xrange(body.N):
				plt.plot(Mesh.x[body.neighbor[k]%Mesh.Nx], Mesh.y[body.neighbor[k]/Mesh.Nx],\
						'ro', markersize=4)
			plt.title('MESH: '+str(Mesh.Nx)+'x'+str(Mesh.Ny)+' / IB: '+str(body.N))
			plt.savefig(Case.path+'/images/mesh_'+str(Mesh.Nx)+'_'+str(Mesh.Ny)+'_'+str(body.N)+'.png')
		else:
			plt.title('MESH: '+str(Mesh.Nx)+'x'+str(Mesh.Ny))
			plt.savefig(Case.path+'/images/mesh_'+str(Mesh.Nx)+'_'+str(Mesh.Ny)+'.png')
		if is_show:
			plt.show()
		plt.clf()
		plt.close()
