# source: $pyIBM/src/timeInfo.py
# Olivier Mesnard (mesnardo@gwu.edu)
# BarbaGroup (lorenabarba.com)

import time


def start():
	'''Starts the timer.'''
	return time.time()

def stop(tic, info=None):
	'''Prints the time spent with a given info.'''
	print '{'+info+'} Execution time: %.3f s' % float(time.time()-tic)
