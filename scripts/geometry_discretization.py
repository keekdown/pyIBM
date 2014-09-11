# file: $pyIBM/scripts/geometry_discretization.py
# author: Olivier Mesnard (mesnardo@gwu.edu)
# BarbaGroup (lorenabarba.com)


import os
import sys
import argparse

import math
import numpy as np
from matplotlib import pyplot as plt


def read_inputs():
	parser = argparse.ArgumentParser(description='Geometry discretization')

	parser.add_argument('-f', '--file', dest='file', type=str,
						help='path of the coordinates file')

	parser.add_argument('-n', dest='n', type=int,
						help='number of points to discretize')

	parser.add_argument('-d', '--delta', dest='delta', type=float,
						help='spatial distance between two consecutive points')

	parser.add_argument('-s', '--scale', dest='scale', type=float,
						help='scale the geometry')

	return parser.parse_args()


class Geometry:
	"""Definition of the geometry."""
	def __init__(self, file_name):
		"""Initializes the geometry by reading the coordinates file.
		
		Arguments
		---------
		file_name -- path of the coordinates file.
		"""
		self.read(file_name)

	def read(self, file_name):
		"""Stores the coordinates of the geometry into arrays.
		
		Arguments
		--------
		file_name -- path of the coordinates file.
		"""
		with open(file_name, 'r') as infile:
			self.x, self.y = np.loadtxt(infile, dtype=float, 
						  	  			delimiter='\t', unpack=True)

	def get_distance(self, x_start, y_start, x_end, y_end):
		"""Returns the distance between two points.
		
		Arguments
		---------
		x_start, y_end -- coordinates of the first point.
		x_end, y_end -- coordinates of the second point.
		"""
		return math.sqrt((x_end-x_start)**2 + (y_end-y_start)**2)

	def get_perimeter(self):
		"""Returns the perimeter of the geometry."""
		x, y = np.append(self.x, self.x[0]), np.append(self.y, self.y[0])
		return np.sum(np.sqrt((x[1:]-x[0:-1])**2+(y[1:]-y[0:-1])**2))

	def get_center_mass(self):
		"""Finds the center of mass."""
		self.x_cm, self.y_cm = self.x.sum()/self.x.size, self.y.sum()/self.y.size

	def rotation(self, aoa=0.0, x_rot=0.0, y_rot=0.0):
		"""Rotates the geometry.
		
		Arguments
		---------
		aoa -- angle of rotation in degrees (default 0.0).
		x_rot, y_rot -- location of the center of rotation (default 0.0, 0.0).
		"""
		aoa *= math.pi/180.
		x_tmp = x_rot + (self.x-x_rot)*math.cos(aoa) - (self.y-y_rot)*math.sin(aoa)
		y_tmp = y_rot + (self.y-y_rot)*math.cos(aoa) - (self.y-y_rot)*math.sin(aoa)
		self.x, self.y = x_tmp, y_tmp
		self.get_center_mass()

	def translation(self, x_trans=0.0, y_trans=0.0):
		"""Translates the geometry.
		
		Arguments
		---------
		x_trans, y_trans -- displacement in the x- and y- directions (default 0.0, 0.0).
		"""
		self.x += x_trans
		self.y += y_trans
		self.get_center_mass()

	def scale(self, ratio=1.0):
		"""Scales the geometry.
		
		Arguments
		---------
		ratio -- scaling ratio (default 1.0).
		"""
		self.x = self.x_cm + ratio*(self.x - self.x_cm)
		self.y = self.y_cm + ratio*(self.y - self.y_cm)
	
	def interpolation(self, x_start, y_start, x_end, y_end, lc):
		"""Computes the coordinates of a point 
		by interpolation between two given points given a distance.
		
		Arguments
		---------
		x_start, y_start -- coordinates of the starting point.
		x_end, y_end -- coordinates of the ending point.
		lc -- length between the starting point and the interpolated one.

		Returns
		-------
		x_target, y_target -- coordinates of the interpolated point.
		"""
		length = self.get_distance(x_start, y_start, x_end, y_end)
		x_target = x_start + lc/length*(x_end-x_start)
		y_target = y_start + lc/length*(y_end-y_start)
		return x_target, y_target

	def projection(self, x_start, y_start, x_tmp, y_tmp, x_end, y_end, lc):
		"""Computes the coordinates of a point
		by projection onto the segment [(x_tmp, y_tmp), (x_end, y_end)]
		such that the distance between (x_start, y_start) and the new point
		is the given cahracteristic length.

		Arguments
		---------
		x_start, y_start -- coordinates of the starting point.
		x_tmp, y_tmp -- coordinates of the intermediate point.
		x_end, y_end -- coordinates of the ending point.
		lc -- characteristic length.

		Returns
		-------
		x_target, y_target -- coordinates of the projected point.
		"""
		# coefficients of the second-order polynomial
		a = (x_end-x_tmp)**2 + (y_end-y_tmp)**2
		b = 2.0*( (x_end-x_tmp)*(y_tmp*(x_start-x_end) + y_end*(x_tmp-x_start)) 
				- y_start*(y_end-y_tmp)**2 )
		c = (y_start**2-lc**2)*(y_end-y_tmp)**2 \
			+ (y_tmp*(x_start-x_end) + y_end*(x_tmp-x_start))**2
		# solve the second-order polynomial: ay^2 + by + c = 0
		y = np.roots([a, b, c])
		# test if the point belongs to the segment
		test = (y_tmp <= y[0] <= y_end or y_end <= y[0] <= y_tmp)
		y_target = (y[0] if test else y[1])
		x_target = x_tmp + (x_end-x_tmp)/(y_end-y_tmp)*(y_target-y_tmp)
		return x_target, y_target

	def discretize(self, N=None, lc=None):
		"""Discretizes the geometry 
		given a characteristic length or a nmumber of points
		
		Arguments
		---------
		N -- number of points (default None).
		lc -- characteristic length (default None).
		"""
		# calculate either the numebr of points or characteristic length
		if N and not lc:
			lc = self.get_perimeter()/N
		elif lc and not N:
			N = int(self.get_perimeter()/lc)

		# exit function if same discretization
		if N == self.x.size:
			return

		# copy coordinates and initialize new ones
		x_old, y_old = np.append(self.x, self.x[0]), np.append(self.y, self.y[0])
		x_new, y_new = np.empty(N, dtype=float), np.empty(N, dtype=float)
		# first element
		x_new[0], y_new[0] = x_old[0], y_old[0]

		I = 0
		tol = 1.0E-06    # tolerance for interpolation
		for i in xrange(N-1):
			x_start, y_start = x_new[i], y_new[i]
			x_end, y_end = x_old[I+1], y_old[I+1]
			distance = self.get_distance(x_start, y_start, x_end, y_end)
			if lc-distance <= tol:
				# interpolation method
				x_new[i+1], y_new[i+1] = self.interpolation(x_start, y_start, x_end, y_end, lc)
			else:
				# projection method
				while I < x_old.size-2 and lc-distance > tol:
					I += 1
					x_tmp, y_tmp = x_end, y_end
					x_end, y_end = x_old[I+1], y_old[I+1]
					distance = self.get_distance(x_start, y_start, x_end, y_end)
				x_new[i+1], y_new[i+1] = self.projection(x_start, y_start, 
														 x_tmp, y_tmp, 
														 x_end, y_end, 
														 lc)
		# stores the new discretization
		self.x, self.y = x_new.copy(), y_new.copy()

	def plot(self):
		"""Plots the geometry."""
		size = 6
		ratio = (self.y.max() - self.y.min()) / (self.x.max() - self.x.min())
		plt.figure(figsize=(size, size*ratio))
		plt.grid(True)
		plt.xlabel(r'$x$', fontsize=18)
		plt.ylabel(r'$y$', fontsize=18)
		plt.plot(np.append(self.x, self.x[0]), 
				 np.append(self.y, self.y[0]),
				 color='b', ls='-', lw=2, marker='o', markersize=6)
		plt.xlim(self.x.min()-1, self.x.max()+1)
		plt.ylim(self.y.min()-1, self.y.max()+1)
		plt.show()

def main():
	# parse the command-line
	args = read_inputs()
	
	# read the coordinates file
	file_name = './scripts/cylinder_test.bdy'
	body = Geometry(file_name)
	body.plot()
	body.translation(x_trans=1.0, y_trans=1.0)
	body.scale(ratio=2.0)
	body.plot()
	body.discretize(N=10)
	body.plot()


if __name__ == '__main__':
	main()
