# source: $pyIMB/src/operations.py
# Olivier Mesnard (mesnardo@gwu.edu)
# BarbaGroup (lorenabarba.com)

def grad(u, direction):
	'''Computes the gradient of a variable in a a given direction.'''
	if direction == 'x':
		return u.gradientx.mat.dot(u.field) + u.gradientx.bc_vect
	elif direction == 'y':
		return u.gradienty.mat.dot(u.field) + u.gradienty.bc_vect

def lap(u):
	'''Computes the laplacian of a variable.'''
	return u.laplacian.mat.dot(u.field) + u.laplacian.bc_vect
