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
	"""A direction contains one or several subdomains."""
	def __init__(self, info_subdomain):
		"""Gets the number of cells or the stretching ratio.
		
		Arguments
		---------
		info_subdomain -- info about the subdomain.
		"""
		self.end = info_subdomain['end']
		if 'cells' in info_subdomain:
			self.N = info_subdomain['cells']
			self.is_uniform = True
		if 'stretchRatio' in info_subdomain:
			self.gamma = info_subdomain['stretchRatio']
			self.is_uniform = False


class Direction:
	"""A mesh contains one or several directions."""
	def __init__(self, info_direction):
		"""Gets info about a direction and creates subdomains.
		
		Arguments
		---------
		info_direction -- info about the direction.
		"""
		self.name = info_direction['direction']
		self.start = info_direction['start']
		self.end = info_direction['subdomains'][-1]['end']
		
		self.subdomain = []
		for info in info_direction['subdomains']:
			self.subdomain.append(Subdomain(info))

	def generate(self):
		"""Generates the mesh in one direction (coordinates and grid spacing)."""
		for idx, sd in enumerate(self.subdomain):
			if sd.is_uniform:
				idx_uniform = idx
				if idx != 0:
					start = self.subdomain[idx-1].end
				else:
					start = self.start
				h = (sd.end-start)/(sd.N-1)
				coord = list(start + h*np.arange(sd.N))
		
		for idx, sd in enumerate(self.subdomain):
			if not sd.is_uniform:
				if idx < idx_uniform:
					while True:
						x = coord[0] + sd.gamma*(coord[0]-coord[1])
						if x-self.start > coord[1]-coord[0]:
							coord.insert(0,x)
						else:
							coord.insert(0,self.start)
							break
				if idx > idx_uniform:
					while True:
						x = coord[-1] + sd.gamma*(coord[-1]-coord[-2])
						if self.end-x > coord[-1]-coord[-2]:
							coord.append(x)
						else:
							coord.append(self.end)
							break
		self.N = len(coord)
		self.coord = np.array(coord)
		self.delta = np.empty(self.N, dtype=float)
		self.delta[0:self.N-1] = self.coord[1:self.N] - self.coord[0:self.N-1]
		self.delta[self.N-1] = self.delta[self.N-2]


class Mesh:
	"""Class to generate a Cartesian grid mesh."""
	def __init__(self):
		"""Parses the file _infoMesh.yaml and generates the mesh."""
		Mesh.is_body = False
		
		infile = open(Case.path+'/_infoMesh.yaml', 'r')
		info_mesh = yaml.load(infile)
		infile.close()
		
		Mesh.direction = []
		for info in info_mesh:
			if 'direction' in info:
				Mesh.direction.append(Direction(info))
			elif 'body' in info:
				Mesh.is_body = True
		
		self.generate()
		
		Mesh.xmin, Mesh.xmax = Mesh.direction[0].start, Mesh.direction[0].end
		Mesh.ymin, Mesh.ymax = Mesh.direction[1].start, Mesh.direction[1].end
		Mesh.Nx, Mesh.Ny = Mesh.direction[0].N, Mesh.direction[1].N
		
		self.write()

		print '\n'
		print '-> Number of points in the x-direction: ', Mesh.Nx
		print '-> Number of points in the y-direction: ', Mesh.Ny
		print '\n'

	def generate(self):
		"""Generates the mesh."""
		for direction in Mesh.direction:
			direction.generate()
		
		Mesh.x, Mesh.y = Mesh.direction[0].coord,Mesh.direction[1].coord
		Mesh.dx, Mesh.dy = Mesh.direction[0].delta,Mesh.direction[1].delta

	def write(self):
		"""Writes the mesh into a data file."""
		with open(Case.path+'/mesh.dat', 'w') as file_name:
			np.savetxt(file_name, np.c_[Mesh.x, Mesh.y, Mesh.dx, Mesh.dy], 
					   fmt='%.6f', delimiter='\t', 
					   header='Mesh (%d by %d): x, y, dx, dy' 
					   % (Mesh.Nx, Mesh.Ny))

	def read(self):
		"""Reads the mesh from a data file."""
		with open(Case.path+'/mesh.dat', 'r') as file_name:
			Mesh.x, Mesh.y, Mesh.dx, Mesh.dy = np.loadtxt(file_name, 
														  dtype=float, 
														  delimiter='\t',
														  unpack=True)
			
	def plot(self,body=None, is_show=False):
		"""Plots the Cartesian mesh grid.
		
		Arguments
		---------
		body -- Body object immersed (default None).
		is_show -- Boolean to display the mesh on the screen (default False).
		"""
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
				plt.plot(Mesh.x[body.neighbor[k]%Mesh.Nx], Mesh.y[body.neighbor[k]/Mesh.Nx],
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
