# source: $pyIMB/src/operations.py
# Olivier Mesnard (mesnardo@gwu.edu)
# BarbaGroup (lorenabarba.com)

def grad(u,direction):
		if (direction=='x'):
			return u.gradient_x.mat.dot(u.field)+u.gradient_x.bc_vect
		elif (direction=='y'):
			return u.gradient_y.mat.dot(u.field)+u.gradient_y.bc_vect

def lap(u):
	return u.laplacian.mat.dot(u.field)+u.laplacian.bc_vect
