# file: $pyIMB/src/operations.py
# author: Olivier Mesnard (mesnardo@gwu.edu)
# BarbaGroup (lorenabarba.com)


def grad(var, direction):
	"""Returns the gradient of a variable in a given direction.
	
	Arguments
	---------
	var -- field variable (Variable object) which gradient is computed.
	direction -- direction to compute the gradient.
	"""
	return ( getattr(var, 'gradient'+direction).mat.dot(var.field)
		   + getattr(var, 'gradient'+direction).bc_vect )

def lap(var):
	"""Returns the Laplacian of a variable.
	
	Arguments
	---------
	var -- field variable (Variable object) which Laplacian is computed.
	"""
	return var.laplacian.mat.dot(var.field) + var.laplacian.bc_vect
