import numpy as np
from math import *
import sys

def main(arg):

	xmin,xmax = -5.0,10.0
	Nx,Ny = 120,120

	R = 0.5
	xc,yc = 0.0,0.0
	Nk = int(2*pi*R*(Nx-1)/(xmax-xmin))
	Xk = np.empty(Nk,dtype=np.float64)
	Yk = np.empty(Nk,dtype=np.float64)
	for k in range(Nk):
		Xk[k] = xc + R*cos(2*pi*k/Nk)
		Yk[k] = yc + R*sin(2*pi*k/Nk)
	outFile = open('./cylinder.bdy','w')
	outFile.write(str(Nk)+'\n')
	for k in range(Nk):
		outFile.write(str(Xk[k])+'\t'+str(Yk[k])+'\n')
	outFile.close()

if (__name__ == '__main__'):
	main(sys.argv)
