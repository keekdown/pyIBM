# file: $pyIBM/src/mesh.py
# author: Olivier Mesnard (mesnardo@gwu.edu)
# BarbaGroup (lorenabarba.com)


import os

import numpy as np
import yaml
from matplotlib import pyplot as plt

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
		
		self.subdomains = []
		for info in info_direction['subdomains']:
			self.subdomains.append(Subdomain(info))

	def generate(self):
		"""Generates the mesh in one direction (coordinates and grid spacing)."""
		for i, subdomain in enumerate(self.subdomains):
			if subdomain.is_uniform:
				i_uniform = i
				if i != 0:
					start = self.subdomains[i-1].end
				else:
					start = self.start
				h = (subdomain.end - start) / (subdomain.N - 1)
				coord = list(start + h*np.arange(subdomain.N))
		
		for i, subdomain in enumerate(self.subdomains):
			if not subdomain.is_uniform:
				if i < i_uniform:
					while True:
						x = coord[0] + subdomain.gamma*(coord[0]-coord[1])
						if x-self.start > coord[1]-coord[0]:
							coord.insert(0, x)
						else:
							coord.insert(0, self.start)
							break
				if i > i_uniform:
					while True:
						x = coord[-1] + subdomain.gamma*(coord[-1]-coord[-2])
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
		
		with open(Case.path+'/_infoMesh.yaml', 'r') as infile:
			info_mesh = yaml.load(infile)
		
		Mesh.directions = {}
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
		
		# write the mesh
		self.write()

		print ('\n')
		print ('-> Number of points in the x-direction: ', Mesh.Nx)
		print ('-> Number of points in the y-direction: ', Mesh.Ny)
		print ('\n')

	def generate(self):
		"""Generates the mesh."""
		for direction in Mesh.directions.values():
			# generate the mesh in each direction
			direction.generate()
	
		# store the mesh coordinates and the cell widths
		Mesh.x, Mesh.y = Mesh.directions['x'].coord, Mesh.directions['y'].coord
		Mesh.dx, Mesh.dy = Mesh.directions['x'].delta, Mesh.directions['y'].delta

	def write(self):
		"""Writes the mesh into a data file."""
		with open(Case.path+'/mesh.dat', 'wb') as outfile:
			np.savetxt(outfile, np.c_[Mesh.x, Mesh.y, Mesh.dx, Mesh.dy], 
					   fmt='%.6f', delimiter='\t', 
					   header='Mesh (%d by %d): x, y, dx, dy' 
					   % (Mesh.Nx, Mesh.Ny))

	def read(self):
		"""Reads the mesh from a data file."""
		with open(Case.path+'/mesh.dat', 'r') as infile:
			Mesh.x, Mesh.y, Mesh.dx, Mesh.dy = np.loadtxt(infile, 
														  dtype=float, 
														  delimiter='\t',
														  unpack=True)
