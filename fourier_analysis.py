#!/usr/bin/env python
import os
from sys import argv
from numpy import *
from matplotlib.pyplot import *
def start_calc(numfiles,Ntot,iso=True):
	hist_file = loadtxt('history.dat')	
	hist = zeros((numfiles,Ntot+1))
	fld = fargo(0)
	sig0 = fld.mrho[:,0];
	hist[0,0] = hist_file[0,0]
	hist[0,1] = 1
	for n in range(1,numfiles):
		print 'Working on file %d out of %d, %.2f%% done' % (n,Ntot,100*float(n)/float(Ntot)) 
		fld=fargo(n,iso)
		hist[n,0] = hist_file[n,0]
		r = fld.r
		norm = trapz(fld.mrho[:,0]*r,x=r)
		val = abs(trapz((fld.mrho[:,0] - sig0)*r,x=r))
		hist[n,1] = norm/val
		for m in range(1,Ntot):
			val = abs(trapz(fld.mrho[:,m]*r,x=r))
			hist[n,m+1] = val/norm
	return hist	


# First argument is directory name
# second argument is number of files to load
# third argument is number modes to log
print 'Importing fargo python module'
execfile('utils/quickreader.py')
dirname = argv[1]
print 'Changing to %s directory' % dirname
os.chdir(dirname)
hist = start_calc(int(argv[2]),int(argv[3]))
print 'Writing results to modehistory.dat'
lines = '\n'.join(['\t'.join([str(x) for x in line]) for line in hist])
with open('modehistory.dat','w') as f:
	f.write(lines)
print 'Finished.'

