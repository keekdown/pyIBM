# file: $pyIBM/src/body.py
# author: Olivier Mesnard (mesnardo@gwu.edu)
# BarbaGroup (lorenabarba.com)
		

from math import sin, cos, pi
import numpy as np
import yaml

from mesh import Mesh
from case import Case


class Body:
	"""Creates an immersed boundary."""
	def __init__(self):
		"""Parses the file _infoBody.yaml 
		and generates the immersed boundary.
		"""
		# parse the file
		with open(Case.path+'/_infoBody.yaml', 'r') as infile:
			info = yaml.load(infile)['IB']
		self.name = info['name']
		self.coord_file = info['coordinates']
		self.is_moving = info['moving']
		# generate the body
		self.generate()

	def generate(self):
		"""Generates the immersed boundary, finds neighbors
		and initializes Lagrangian variables.
		"""
		# reads the coordinate file
		with open(Case.path+'/'+self.coord_file, 'r') as file_name:
			self.x, self.y = np.loadtxt(file_name, dtype=float, 
										delimiter='\t', unpack=True)
		
		self.N = len(self.x)    # number of body points
		print '\n-> Number of points on the body: ', self.N, '\n'
		
		# computes the length vector
		self.dx = np.empty(self.N, dtype=float)
		self.dy = np.empty(self.N, dtype=float)
		self.dx[1:self.N] = self.x[1:self.N] - self.x[0:self.N-1]
		self.dy[1:self.N] = self.y[1:self.N] - self.y[0:self.N-1]
		self.dx[0] = self.x[0] - self.x[-1]
		self.dy[0] = self.y[0] - self.y[-1]
	
		# find the neighbors
		self.get_neighbors()
		
		# initializes other Lagrangian variables
		self.u = np.empty(self.N, dtype=float)
		self.v = np.empty(self.N, dtype=float)
		self.fx = np.empty(self.N, dtype=float)
		self.fy = np.empty(self.N, dtype=float)
		self.ud = np.zeros(self.N, dtype=float)
		self.vd = np.zeros(self.N, dtype=float)
		self.x0, self.y0 = np.copy(self.x), np.copy(self.y)

	def get_neighbors(self):
		"""Finds the closest Eulerian point on the mesh grid."""
		self.neighbor = np.empty(self.N, dtype=int)
		for k, (x, y) in enumerate(zip(self.x, self.y)):
			for i in xrange(Mesh.Nx-1):
				if Mesh.x[i] <= x < Mesh.x[i+1]:
					I = i
					break
			for j in xrange(Mesh.Ny-1):
				if Mesh.y[j] <= y < Mesh.y[j+1]:
					J = j
					break
			self.neighbor[k] = J*Mesh.Nx + I

	def kinematics(self):
		"""Defines the kinematics of a moving immersed boundary."""
		A = 1.0
		Nf = 100
		f = 1./(Nf*Solver.dt)
		self.x[:] = self.x0[:] + 0.0
		self.x[:] = self.y0[:] + A*sin(2*pi*f*Solver.ite*Solver.dt)
		self.ud[:] = 0.0
		self.vd[:] = 2*pi*f*A*cos(2*pi*f*Solver.ite*Solver.dt)
		self.getNeighbors()
