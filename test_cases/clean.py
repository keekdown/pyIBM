import sys
import os

def main(arg):

	print '[info]: cleaning solution folders in ',arg[1]
	os.system('rm -rf '+arg[1]+'/0*')
	os.system('rm -rf '+arg[1]+'/1*')
	os.system('rm -rf '+arg[1]+'/2*')
	os.system('rm -rf '+arg[1]+'/3*')
	os.system('rm -rf '+arg[1]+'/4*')
	os.system('rm -rf '+arg[1]+'/5*')
	os.system('rm -rf '+arg[1]+'/6*')
	os.system('rm -rf '+arg[1]+'/7*')
	os.system('rm -rf '+arg[1]+'/8*')
	os.system('rm -rf '+arg[1]+'/9*')

	print '[info]: cleaning pictures in ',arg[1]
	os.system('rm -rf '+arg[1]+'/images')
	os.system('rm -f '+arg[1]+'/mesh*.dat')
	os.system('rm -f '+arg[1]+'/*.png')

if (__name__=='__main__'):
	print '\n\t----- START - CLEANING -----\n'
	main(sys.argv)
	print '\n\t----- END - CLEANING -----\n'
