#!/usr/bin/env python
################################################################################
#Voronoi Tessellation with Fortune's sweep line algorithm
#Author: Teng Liu
#Email: lewtonstein@gmail.com
################################################################################
import numpy as np
from astropy.io import fits
import json
from SweeplineVT import Voronoi
import time,sys,warnings,os,getopt

def main():
	#	Notice the input coordinate is taken as in Python (image: 0~x-1,0~y-1; event: -0.5~x-0.5,-0.5~y-0.5)
	#	rather than the ds9 style (image: 1~y,1~x; event: 0.5~y+0.5,0.5~x+0.5).
	#-P,--calpvd		Calculate Pixelated Voronoi Diagram, only in case of image input.
	#-T,--calTriangle	Divide the space into smallest triangles.
	#-S,--caldst		Calculate Density image, only in case of image input.
	#-D,--calDelaunay	Make Delaunay diagram.
	#--resolution N		Set position resolution to the Nth decimal place.
	#			Only in case of events input, where high resolution can make problems.
	#--makeimage		Make an image from events list
	#			The image size rests with the max coordinate.
	#			Rescale the input positions by yourself.
	#			Use log scale if necessary.
	#			No need to reset border with --border, unless you want to enlarge the image.
#Caveat
#	PVD (Pixelated Voronoi Diagram) is not uniquely defined. There is uncertainty in some pixels.
#
#	In the SAOds9 reg file:
#		X_ds9, Y_ds9 = Y_py+1, X_py+1
#		the image border (range of the Voronoi diagram) is -0.5~size-0.5 for X_py, Y_py, and 0.5~size+0.5 for X_ds9, Y_ds9.
#	The way to reload the Voronoi diagram reg file given by this program is:
#		np.fromregex('xxx.reg',"line\(([0-9.]*),([0-9.]*),([0-9.]*),([0-9.]*)\) # tag={([0-9.]*),([0-9.]*),([0-9.]*),([0-9.]*)}\\n",np.float32)

	__doc__="""
Mode 1:
	slvt.py File [options]
Mode 2:
	slvt.py Number --makeCVT --border x1,x2,y1,y2

The input File can be:
	an image
		If the file has a ".fits" suffix.
	OR a list of coordinates
		It should be an ASCII file containing the point coordinates in its first two columns.

OPTIONS
	-A,--calArea		Calculate cell Areas
	-C,--calCentroid	Calculate cell centroid
	-M,--makeCVT	Make Centroidal Voronoi Tessellation
	--border x1,x2,y1,y2	Set image border in case of events input.
                                e.g. -0.5,1023.5,-0.5,1023.5 --> a 1024x1024 image to view in ds9
	--rmedgepoint		Remove edge points in Area map
	-h/--help		Help

NOTE
	When --calCentroid or --makeCVT, you might want to specify --border explicitly, or else extra padding will be added just for better looking.
		"""
	def usage():
		print(__doc__)
		exit()
	Options={}
	S_opt='dDAPTSMhs'
	L_opt=['calpvd','calArea','calCentroid','caldst','calDelaunay','calTriangle','rmedgepoint','makeCVT','CleanFilledBorder','border=','resolution=','accuracy=','makeimage','help','silent','noautoscale']
	opts,args=getopt.getopt(sys.argv[1:],S_opt,L_opt)
	if len(args)>0:
		for arg in args:
			if os.path.isfile(arg):
				InputFile=arg
				sys.argv.remove(arg)
				break
			elif arg.isdigit():
				CVTNumber=int(arg)
				sys.argv.remove(arg)
				break
		opts,args=getopt.getopt(sys.argv[1:],S_opt,L_opt)
	for opt,arg in opts:
		if opt == '--calpvd' or opt == '-P':
			Options['calpvd']=True
		elif opt == '--calDelaunay' or opt == '-D':
			Options['calDelaunay']=True
		elif opt == '--calTriangle' or opt == '-T':
			Options['calTriangle']=True
		elif opt == '--calArea' or opt == '-A':
			Options['calArea']=True
		elif opt == '--calCentroid' or opt == '-C':
			Options['calCentroid']=True
		elif opt == '--makeCVT' or opt == '-M':
			Options['calCentroid']=True
			Options['makeCVT']=True
		elif opt == '--caldst' or opt == '-S':
			Options['caldst']=True
		elif opt == '--silent' or opt == '-s':
			Options['Silent']=True
			warnings.simplefilter('ignore')
		elif opt == '--noautoscale':
			Options['autoscale']=False
		elif opt == '--makeimage':
			Options['MakeIntImage']=True
		elif opt == '--rmedgepoint':
			Options['RemoveEdgePoint']=True
		elif opt == '--CleanFilledBorder':
			Options['CleanFilledBorder']=True
		elif opt == '--border':
			try:
				arg = arg.split(',')
				assert len(arg)==4
				border={}
				if arg[0]!='': border['xlow']=float(arg[0])
				if arg[1]!='': border['xhigh']=float(arg[1])
				if arg[2]!='': border['ylow']=float(arg[2])
				if arg[3]!='': border['yhigh']=float(arg[3])
			except:
				print('ERROR: --border',arg)
				exit()
			else:
				Options['border'] = border
		elif opt == '--resolution':
			try:
				n = int(arg)
				assert n>=0
			except:
				print('ERROR: --resolution',arg)
				exit()
			else:
				if n>=Voronoi.SLVround: raise RuntimeError('--resolution Do you really want such a high resolution?')
				Options['Resolution'] = n
		elif opt == '--accuracy':
			try:
				n = int(arg)
				assert n>0
			except:
				print('ERROR: --accuracy',arg)
				exit()
			else:
				Voronoi.SLVround = n
		elif opt == "-h" or opt == "--help":
			usage()
		elif opt == '-d':
			Voronoi.debug = True
	if len(args)>0: sys.exit("I don't understand"+str(args))
	if 'InputFile' in locals():
		Options['FileName']=InputFile.rsplit('.',1)[0]
		if len(InputFile.rsplit('.',1))>1 and InputFile.rsplit('.',1)[1] == 'fits':
			data,hdr=fits.getdata(InputFile,header=True)
			Options['Hdr']=hdr
			vor=Voronoi(image=data,**Options)
		else: #a file which store the coordinates of points in the first two columns
			data=np.loadtxt(InputFile)
			vor=Voronoi(events=data,**Options)
			CVTNumber=len(data)
	elif 'CVTNumber' in locals():
		if Options.get('makeCVT',False) and Options.get('border',False):
			data=np.random.random(size=(CVTNumber,2))*[border['xhigh']-border['xlow'],border['yhigh']-border['ylow']]+[border['xlow'],border['ylow']]
			vor=Voronoi(events=data,**Options)
			assert vor.RightLimit+0.5==border['xhigh']-border['xlow'] \
			and vor.TopLimit+0.5==border['yhigh']-border['ylow']
		else: sys.exit('A number is accepted only in makeCVT mode with given border')
	else:
		sys.exit("Please input an image, or a list of points, or a number!\n")

	if not Options.get('makeCVT',False):
		vor.saveresults()
	else:
		ctd=np.array(list(vor.Wmap.values()))
		xoff,yoff=vor.OffSetX,vor.OffSetY
		border={'xlow':-0.5,'ylow':-0.5,'xhigh':vor.RightLimit,'yhigh':vor.TopLimit}
		scale=vor.scale
		MaxIteration=1000
		MaxSep=0.0001
		for n in range(MaxIteration):
			Options['FileName']='iter'+str(n)
			Options['border']=border
			if np.any(np.isnan(ctd)):
				print(ctd)
				exit()
			vor=Voronoi(events=ctd,**Options)
			ctd=np.array(list(vor.Wmap.values()))
			#if n%10==0: vor.saveresults()
			#print('Max',np.max(np.sqrt(np.sum(vor.ctdvector**2,axis=1))))
			if np.max(np.sqrt(np.sum(vor.ctdvector**2,axis=1)))<MaxSep: break
		if n==MaxIteration-1: warnings.warn(f'Limited to {MaxIteration} iteration. Not yet optimized')
		Options['FileName']='CVT'+str(CVTNumber)
		vor=Voronoi(events=ctd,**Options)
		vor.OffSetX=xoff
		vor.OffSetY=yoff
		vor.saveresults()

#import profile
if __name__ == '__main__':
	#profile.run('main()')
	main()
