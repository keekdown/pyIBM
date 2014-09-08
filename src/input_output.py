# file: $pyIBM/src/input_output.py
# author: Oliver Mesnard (mesnardo@gwu.edu)
# BarbaGroup (lorenabarba.com)


import time
import argparse


def read_inputs():
	# list of command-line arguments
	parser = argparse.ArgumentParser(description='Solves the Navier-Stokes equations')
	parser.add_argument('-p', '--path', dest='path', type=str,
						help='path of the case folder')
	parser.add_argument('--mesh', dest='mesh', action='store_true',
						help='generates the computational domain without solving the Navier-Stokes equations')

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
	print '{%s} Execution time: %.3f s' % (info, float(time.time()-tic))
