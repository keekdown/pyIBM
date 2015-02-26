# file: $pyIBM/src/input_output.py
# author: Oliver Mesnard (mesnardo@gwu.edu)
# BarbaGroup (lorenabarba.com)


import time
import argparse


def read_inputs():
	"""Parser of pyIBM."""
	# create the parser
	parser = argparse.ArgumentParser(description='pyIBM command-line arguments')
	
	# fill the parser
	
	parser.add_argument('-p', '--path', dest='path', type=str,
						help='path of the case folder')

	parser.add_argument('-f', '--file', dest='file', type=str,
						help='path of the file')

	parser.add_argument('--mesh', dest='mesh', action='store_true',
						help=('generates the computational domain ' + 
							 'without solving the Navier-Stokes equations'))
	
	parser.add_argument('-v', '--variable', dest='variable', type=str, nargs='+',
						default=['pressure', 'velocity', 'vorticity'],
					    help='list of flow variables to plot')

	parser.add_argument('-z', '--zoom', dest='zoom', type=float, nargs='+',
						help='sets the limits of the figures')

	parser.add_argument('-t', '--time', type=int, nargs='+',
						help=('time-level(s) to plot ' + 
							  '(either a list [min, max, every] or a specific time-level)'))

	# parse the command-line
	return parser.parse_args()


def timer_start():
	"""Returns the time."""
	return time.time()


def timer_stop(tic, info=None):
	"""Stops timer and prints the time with a given info.
	
	Arguments
	---------
	tic -- time when the timer starts.
	info -- info to ouput (default None).
	"""
	print ('{%s} Execution time: %.3f s' % (info, float(time.time()-tic)))
