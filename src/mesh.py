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
		infile = open(Case.path+'/_infoMesh.yaml', 'r')
		info_mesh = yaml.load(infile)
		infile.close()
		
		Mesh.directions = {}
		Mesh.is_body = False
		for info in info_mesh:
			if 'direction' in info:
				Mesh.directions[info['direction']] = Direction(info)
			elif 'body' in info:
				Mesh.is_body = True
		
		# generate the mesh
		self.generate()
		
		# store the boundaries of the domain and the number of points in each direction
		Mesh.xmin, Mesh.xmax = Mesh.directions['x'].start, Mesh.directions['x'].end
		Mesh.ymin, Mesh.ymax = Mesh.directions['y'].start, Mesh.directions['y'].end
		Mesh.Nx, Mesh.Ny = Mesh.directions['x'].N, Mesh.directions['y'].N
		
		self.write()

		print '\n'
		print '-> Number of points in the x-direction: ', Mesh.Nx
		print '-> Number of points in the y-direction: ', Mesh.Ny
		print '\n'

	def generate(self):
		"""Generates the mesh."""
		for direction in Mesh.directions.itervalues():
			# generate the mesh in each direction
			direction.generate()
	
		# store the mesh coordinates and the cell widths
		Mesh.x, Mesh.y = Mesh.directions['x'].coord, Mesh.directions['y'].coord
		Mesh.dx, Mesh.dy = Mesh.directions['x'].delta, Mesh.directions['y'].delta

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
			
	def plot(self, body=None, is_show=False):
		"""Plots the Cartesian mesh grid.
		
		Arguments
		---------
		body -- Body object immersed in the domain (default None).
		is_show -- Boolean to display the mesh on the screen (default False).
		"""
		if not os.path.isdir(Case.path+'/images'):
			os.system('mkdir '+Case.path+'/images')
		plt.figure(num=None)
		plt.grid(False)
		plt.xlabel(r'$x$', fontsize=16)
		plt.ylabel(r'$y$', fontsize=16)
		for i in xrange(Mesh.Nx):
			plt.axvline(Mesh.x[i])
		for j in xrange(Mesh.Ny):
			plt.axhline(Mesh.y[j])
		plt.xlim(Mesh.xmin, Mesh.xmax)
		plt.ylim(Mesh.ymin, Mesh.ymax)
		title = 'MESH: %d x %d' % (Mesh.Nx, Mesh.Ny)
		save_name = 'mesh_%d_%d' % (Mesh.Nx, Mesh.Ny)
		if body != None:
			plt.plot(np.append(body.x, body.x[0]), np.append(body.y, body.y[0]),
					 color='k', ls='-', lw=1, marker='o', markersize=4)
			title += ' / IB: %d' % body.N
			save_name += '_%d' % body.N
		plt.title(title, fontsize=16)
		plt.savefig(Case.path+'/images/'+save_name+'.png')
		if is_show:
			plt.show()
		plt.clf()
		plt.close()
