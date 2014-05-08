# source: $pyIBM/src/timeInfo.py
# Olivier Mesnard (mesnardo@gwu.edu)
# BarbaGroup (lorenabarba.com)

import time

'''Routines to get the time.'''

def start():
	return time.time()

def stop(tic,info=None):
	print '{'+info+'} Execution time: ',('%.3f'%float(time.time()-tic)),' s'
