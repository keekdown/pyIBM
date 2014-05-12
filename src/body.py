# source: $pyIBM/src/body.py
# Olivier Mesnard (mesnardo@gwu.edu)
# BarbaGroup (lorenabarba.com)
		
from math import *
import numpy as np
import yaml

from mesh import Mesh

class Body:
	'''Creates an immersed boundary.'''

	def __init__(self, info_body):
		self.info_body = info_body
		infile = open(self.info_body, 'r')
		info = yaml.load(infile)
		infile.close()
		self.name = info['IB']['name']
		self.coord_file = info['IB']['coordinates']
		self.is_moving = info['IB']['moving']
		self.generate()
		self.cl, self.cd = 0., 0.

	def generate(self):
		'''Generate the immersed boundary,
		and initialize Lagrangian variables.
		'''
		infile = open(Mesh.case_path+'/'+self.coord_file, 'r')
		index, i = 0, 0
		for line in infile:
			data = line.split()
			if index == 0:
				self.N = int(data[0])
				self.x = np.empty(self.N, dtype=float)
				self.y = np.empty(self.N, dtype=float)
				self.dx = np.empty(self.N, dtype=float)
				self.dy = np.empty(self.N, dtype=float)
				self.neighbor = np.empty(self.N, dtype=int)
				index += 1
			else:
				self.x[i], self.y[i] = float(data[0]), float(data[1])
				i += 1
		infile.close()
		self.dx[1:self.N] = self.x[1:self.N] - self.x[0:self.N-1]
		self.dy[1:self.N] = self.y[1:self.N] - self.y[0:self.N-1]
		self.dx[0] = self.x[0] - self.x[-1]
		self.dy[0] = self.y[0] - self.y[-1]
		self.get_neighbors()
		print '\n-> Number of points on the body: ', self.N, '\n'
		self.u = np.empty(self.N, dtype=float)
		self.v = np.empty(self.N, dtype=float)
		self.fx = np.empty(self.N, dtype=float)
		self.fy = np.empty(self.N, dtype=float)
		self.ud = np.zeros(self.N, dtype=float)
		self.vd = np.zeros(self.N, dtype=float)
		self.x0, self.y0 = np.copy(self.x), np.copy(self.y)

	def get_neighbors(self):
		'''Find the closest Eulerian points,
		for each Lagrangian point on the immersed boundary.
		'''
		for k in xrange(self.N):
			for i in xrange(Mesh.Nx-1):
				if Mesh.x[i] <= self.x[k] < Mesh.x[i+1]:
					I = i
			for j in xrange(Mesh.Ny-1):
				if Mesh.y[j] <= self.y[k] < Mesh.y[j+1]:
					J = j
			self.neighbor[k] = J*Mesh.Nx+I

	def kinematics(self):
		'''Define the kinematics of a moving immersed boundary.
		This function needs to be adapted to a given problem/application.
		'''
		A = 1.0
		Nf = 100
		f = 1./(Nf*Solver.dt)
		self.x[:] = self.x0[:] + 0.0
		self.x[:] = self.y0[:] + A*sin(2*pi*f*Solver.ite*Solver.dt)
		self.ud[:] = 0.0
		self.vd[:] = 2*pi*f*A*cos(2*pi*f*Solver.ite*Solver.dt)
		self.getNeighbors()
