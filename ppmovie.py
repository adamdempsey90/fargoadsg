#!/usr/bin/env python
import os
from sys import argv
import matplotlib
matplotlib.use('Agg')
def save_images(q,filelist,iso=True,dirname='./',fnamebase='image',imgext='.png'):
	if dirname[-1] != '/':
		dirname += '/'
	totnum = len(filelist)
	for i,j in enumerate(filelist):
		fld = fargo(j,iso)
		fname = dirname + fnamebase + '%03d'%i + imgext
		print 'Saving image %d/%d to %s...\t %.2f%% done' % (i,totnum,fname,float(i)/float(totnum)*100)
		fld.plot(q,output=True,fname=fname)
	


# First argument is directory name
# second argument is quantity to plot, see fargo class plot method for values
# third argument is number of files to load
# fourth argument is total number of dumps
# fifth argument is output directory in first argument directory
# sixth argument is extension, defaults to png if not given
numargs = len(argv)

print 'Importing fargo python module'
execfile('utils/quickreader.py')
dirname = argv[1]
q = argv[2]
numfiles = int(argv[3])
totfiles = int(argv[4])
imgdir = argv[5]
if numargs == 7:
	ext = argv[6]
else:
	ext = '.png'




filelist = range(totfiles)[::totfiles/numfiles]
print filelist
print 'Changing to %s directory' % dirname
os.chdir(dirname)
try:
	print 'Making %s directory' % imgdir
	os.mkdir(imgdir)
except:
	pass

save_images(q,filelist,dirname=imgdir,fnamebase=q,imgext=ext);
print 'Finished.'

