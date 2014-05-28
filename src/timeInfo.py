# source: $pyIBM/src/timeInfo.py
# Olivier Mesnard (mesnardo@gwu.edu)
# BarbaGroup (lorenabarba.com)

import time


def start():
	"""Returns the time."""
	return time.time()

def stop(tic, info=None):
	"""Stops timer and prints the time with a given info.
	
	Arguments
	---------
	tic -- previous time to get time spent.
	info -- information to display (default None).
	"""
	print '{'+info+'} Execution time: %.3f s' % float(time.time()-tic)
