# source: $pyIMB/src/operations.py
# Olivier Mesnard (mesnardo@gwu.edu)
# BarbaGroup (lorenabarba.com)

def grad(u, direction):
	'''Compute the gradient of 'u' in a a given direction.'''

	if direction == 'x':
		return u.gradientx.mat.dot(u.field) + u.gradientx.bc_vect
	elif direction == 'y':
		return u.gradienty.mat.dot(u.field) + u.gradienty.bc_vect

def lap(u):
	'''Compute the laplacian of a vector.'''

	return u.laplacian.mat.dot(u.field) + u.laplacian.bc_vect
