# source: $pyIMB/src/operations.py
# Olivier Mesnard (mesnardo@gwu.edu)
# BarbaGroup (lorenabarba.com)

def grad(var, direction):
	'''Computes the gradient of a variable in a given direction.'''
	return ( getattr(var, 'gradient'+direction).mat.dot(var.field)
		   + getattr(var, 'gradient'+direction).bc_vect )

def lap(var):
	'''Computes the laplacian of a variable.'''
	return var.laplacian.mat.dot(var.field) + var.laplacian.bc_vect
