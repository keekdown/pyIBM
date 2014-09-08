# file : $pyIBM/src/domain.py
# author: Olivier Mesnard (mesnardo@gwu.edu)
# BarbaGroup (lorenabarba.com)


import sys

import numpy as np
from matplotlib import pyplot as plt

from case import Case
from mesh import Mesh
from body import Body


class Domain:
	"""A Domain contains a Mesh and possibly an immersed boundary."""
	def __init__(self, is_show=False):
		"""Creates the domain by generating a mesh
		and possibly an immersed boundary.
		
		Arguments
		---------
		is_show -- Boolean to show the domain (default False).
		"""
		# generate the mesh
		Domain.mesh = Mesh()
		Domain.is_body = Mesh.is_body

		# generate the immersed boundary if one
		Domain.body = None
		if Domain.is_body:
			Domain.body = Body()

		# plot the domain
		self.plot(is_show)

	def plot(self, is_show):
		"""Plots the computational domain.
		
		Arguments
		---------
		is_show -- Boolean to show the domain on the screen.
		"""
		size= 6
		ratio = (Mesh.ymax - Mesh.ymin) / (Mesh.xmax - Mesh.xmin)
		plt.figure(figsize=(size, ratio*size))
		plt.xlabel(r'$x$', fontsize=18)
		plt.ylabel(r'$y$', fontsize=18)
		# plot the mesh
		for x in Mesh.x:
			plt.axvline(x)
		for y in Mesh.y:
			plt.axhline(y)
		plt.xlim(Mesh.xmin, Mesh.xmax)
		plt.ylim(Mesh.ymin, Mesh.ymax)
		title = 'MESH: %d x %d' % (Mesh.Nx, Mesh.Ny)
		save_name = 'mesh_%d_%d' % (Mesh.Nx, Mesh.Ny)
		if Domain.is_body:
			# plot the immersed boundary
			plt.plot(np.append(Domain.body.x, Domain.body.x[0]),
					 np.append(Domain.body.y, Domain.body.y[0]),
					 color='k', ls='-', lw=1.0, marker='o', markersize=4)
			title += ' / IB: %d' % Domain.body.N
			save_name += '_%d' % Domain.body.N
		plt.title(title, fontsize=14)
		plt.savefig(Case.path+'/images/'+save_name+'.png')
		if is_show:
			plt.show()
			sys.exit(0)
		plt.clf()
		plt.close()
