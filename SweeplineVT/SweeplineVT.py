#!/usr/bin/env python
################################################################################
#Voronoi Tessellation with Fortune's sweep line algorithm
#Author: Teng Liu
#Email: lewtonstein@gmail.com
################################################################################
import numpy as np
from matplotlib import pyplot as pl
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from astropy.io import fits
import heapq,json
import itertools
from copy import copy
import time,sys,warnings
from tqdm import tqdm
from numba import jit,boolean,float64,int32
from scipy.interpolate import griddata
def _ltshowwarning(message, category, filename, lineno, file=None, line=None):
	if category==UserWarning: print('\033[33mWARNING\033[0m:',message)
	else:
		msg = warnings.WarningMessage(message, category, filename, lineno, file, line)
		warnings._showwarnmsg_impl(msg)
warnings.showwarning=_ltshowwarning
color = lambda s,ncolor,nfont: "\033["+str(nfont)+";"+str(ncolor)+"m"+s+"\033[0;0m"
#define:
#p<q if py<qy or py=qy and px<qx
#vertical,upward: y++
#horizontal,right: x++

def INT2RGB(RGBint):
    return np.array([(RGBint >> 16) & 255,(RGBint >> 8) & 255,RGBint & 255])/255.
def INT2HEX(RGBint):
    return f'#{(RGBint >> 16)&255:02x}{(RGBint >> 8)&255:02x}{RGBint&255:02x}'
def RGB2INT(rgb):
    return ((rgb[0]<<16) + (rgb[1]<<8) + rgb[2])
def imgRGB2INT(rgb):
    return ((rgb[:,:,0]<<16) + (rgb[:,:,1]<<8) + rgb[:,:,2])

if True:
	from CweeplineVT import CalCxRightOf,CalIntersect
else:
	@jit(float64(float64,float64,float64,float64,int32,float64,float64,float64),nopython=True,cache=True,nogil=True,fastmath=False)
	def CalCxRightOf(u0,u1,v0,v1,self_direct,px,y0,Voronoi_atol):
		if abs(u1-v1)<Voronoi_atol: cx = (u0+v0)/2.  #np.isclose is very slow
		elif abs(u1-y0)<Voronoi_atol: cx = u0
		elif u0==v0==px: cx=np.float64(np.inf)*self_direct #fastmath=True doesn't work with inf
		else:
			insqrt = (u1*v1-u1*y0-v1*y0+y0*y0)*(u0*u0-2*u0*v0+u1*u1-2*u1*v1+v0*v0+v1*v1)
			cx = (self_direct*np.sqrt(insqrt)-u0*v1+u0*y0+u1*v0-v0*y0)/(u1-v1)
			if insqrt<0: print('Error sqrt<0',insqrt,u0,u1,v0,v1,px,y0,cx)
		return cx #has to be sent back to the Python object to set the hist=(y0,cx)
	#END_OF_def CalCxRightOf

	@jit(float64[:](float64,float64,float64,float64,float64,float64),nopython=True,cache=True,nogil=True,fastmath=False)
	def CalIntersect(ux,uy,vx,vy,wx,wy):
		if uy == vy and ux != vx and vy != wy:
			cx = (ux + vx)/2
			if wx==cx:
				cy = (uy+wy)/2 - (ux-vx)**2/8./(wy-uy)
				ytop=wy
			else:
				#solve[(x - u1)^2+(y - u2)^2=r2,(x-v1)^2+(y-u2)^2=r2, (x-w1)^2+(y-w2)^2=r2, x,y]
				#cy = (-ux*wx - vx*wx + wx**2 + wy**2 + ux*vx - uy**2)/2./(wy-uy)
				cy = (wy+uy)/2. +  (ux-wx)*(vx-wx)/2./(wy-uy)
				r = np.sqrt(((ux-wx)**2+(uy-wy)**2)*((wx-vx)**2+(uy-wy)**2))/2./(wy-uy)
				ytop=cy+r
		elif wx==ux and wx!=vx:
			cy = (wy+uy)/2
			if vy==cy:
				cx = (ux+vx)/2 + (uy-wy)**2/(ux-vx)/8.
				r = ((uy-wy)**2+4*(ux-vx)**2)/8./abs(ux-vx)
				ytop=cy+r
			else:
				cx = (ux+vx)/2 + (vy-wy)*(uy-vy)/(ux-vx)/2 # 2(x-(ux+vx)/2) / (vy-wy) = (uy-vy) / (ux-vx) #same angle
				r = np.sqrt(((ux-vx)**2+(uy-vy)**2)*((wy-vy)**2+(ux-vx)**2))/2./abs(ux-vx)
				ytop=cy+r
		else:
			a = (ux-vx)*(vy-wy)-(uy-vy)*(vx-wx)
			if a==0: return np.float64([-np.inf,np.inf,np.inf]) #take it as None !!!!!!!!!!!!!!check
			b1 = (ux-vx)*(ux+vx)+(uy-vy)*(uy+vy)
			b2 = (vx-wx)*(vx+wx)+(vy-wy)*(vy+wy)
			cx = (b1*(vy-wy)-b2*(uy-vy))/2./a
			cy = (b2*(ux-vx)-b1*(vx-wx))/2./a
			r = np.sqrt((cx-ux)**2+(cy-uy)**2)
			ytop=cy+r
		return np.float64([cx,ytop,cy])
	#END_OF_def CalIntersect

class HalfEdgeM(object):
	#RightLimit=0
	#TopLimit=0
	def __init__(self,p0,p1,direct,base,twin=None,summit=None):
		#if p0[1]<p1[1] or p0[1]==p1[1] and p0[0]<p1[0]:
		#	#the 2 points on 2 sides, p0 < p1
		self.p0 = p0
		self.p1 = p1
		self.direct = direct # 0:upward -1:left 1:right
		self.base = base
		#history test, to save time
		if self.direct == 0:
			self.hist = ((self.p0[0]+self.p1[0])/2.,(self.p0[0]+self.p1[0])/2.)
		else:
			self.hist = (-np.inf,np.inf)
		self.twin = twin
		self.summit = summit
		#self.length = None
		#self.PPdist = None
	def __str__(self):
		return 'p0:'+str(self.p0)+' p1:'+str(self.p1)+' dir:'+str(self.direct)+' b:'+str(self.base)+' s:'+str(self.summit)+' hist:'+str(self.hist)

	def isRightOf(self,px,y0,Voronoi_atol):
		#check whether this edge is to the right of the site px,y0 (the intersection x-position cx > px)
		#(px,y0) is a site event, so y0 is the current position of sweep line
		#Each position of the sweep line defines a circle crossing two neighboring sites (self_p0 and self_p1) and touching the sweep line at the top. As y0 grows, the center of the circle (cx,y0-r) traces the edge.
		#(the method from the paper doesn't work, so I use mine)
		#Input the following expression in http://www.wolframalpha.com/
		#solve[(x-u0)^2+(y-u1)^2=r^2,(x-v0)^2+(y-v1)^2=r^2,y0-y=r,x,y,r]
		#to get cx as a function of self_p0, self_p1, y0
		#There are two solutions, indicating the two crossing points of two nerghboring parabolas
		#The first crossing point (the one nearer to the flatter parabola) is what we want.
		if self.p0 == (-1,-1) : return False
		if self.p1 == (-2,-2): return True
		if self.direct == 0: cx = self.hist[1]
		else:
			if self.hist[0] == y0: cx = self.hist[1]
			else:
				if self.p1[1]>self.p0[1]: u,v = self.p1,self.p0
				else: v,u = self.p1,self.p0
				if u[1]>y0: sys.exit('Error u1>y0   %.8f %.8f - %.8f %.8f - %.8f %.8f' % (u[0],u[1],v[0],v[1],px,y0))
				cx=np.float64(CalCxRightOf(u[0],u[1],v[0],v[1],self.direct,px,y0,Voronoi_atol))
				if not np.isinf(cx): self.hist=[y0,cx]
			if Voronoi.debug: print('isRightOf',self.p0,self.p1,self.direct,px,y0,cx)
		return cx>px

	def oldisRightOf(self,px,y0,Voronoi_atol=1e-8):
		#check whether this edge is to the right of the site px,y0 (the intersection x-position cx > px)
		#(px,y0) is a site event, so y0 is the current position of sweep line
		#Each position of the sweep line defines a circle crossing two neighboring sites (self.p0 and self.p1) and touching the sweep line at the top. As y0 grows, the center of the circle (cx,y0-r) traces the edge.
		#(the method from the paper doesn't work, so I use mine)
		#Input the following expression in http://www.wolframalpha.com/
		#solve[(x-u0)^2+(y-u1)^2=r^2,(x-v0)^2+(y-v1)^2=r^2,y0-y=r,x,y,r]
		#to get cx as a function of self.p0, self.p1, y0
		#There are two solutions, indicating the two crossing points of two nerghboring parabolas
		#The first crossing point (the one nearer to the flatter parabola) is what we want.
		if self.p0 == (-1,-1) : return False
		if self.p1 == (-2,-2): return True
		if self.direct == 0: return self.hist[0]>px
		else:
			if self.hist[0] == y0:
				#print('hist',y0,self.hist)
				cx = self.hist[1]
			else:
				if self.p1[1]>self.p0[1]: u,v = self.p1,self.p0
				else: v,u = self.p1,self.p0
				#print(u,v,y0)
				if u[1]>y0: sys.exit('Error u1>y0   %.8f %.8f - %.8f %.8f - %.8f %.8f' % (u[0],u[1],v[0],v[1],px,y0))
				if round(u[1],Voronoi.SLVround)==round(v[1],Voronoi.SLVround): cx = (u[0]+v[0])/2.
				elif round(u[1],Voronoi.SLVround)==round(y0,Voronoi.SLVround): cx = u[0]
				elif u[0]==v[0]==px: cx=np.inf*self.direct
				else:
					insqrt = (u[1]*v[1]-u[1]*y0-v[1]*y0+y0*y0)*(u[0]*u[0]-2*u[0]*v[0]+u[1]*u[1]-2*u[1]*v[1]+v[0]*v[0]+v[1]*v[1])
					cx = (self.direct*np.sqrt(insqrt)-u[0]*v[1]+u[0]*y0+u[1]*v[0]-v[0]*y0)/(u[1]-v[1])
					if insqrt<0: print('Error sqrt<0 %g %.8f %.8f - %.8f %.8f - %.8f %.8f  %.8f' % (insqrt,u[0],u[1],v[0],v[1],px,y0,cx))
				if not np.isinf(cx): self.hist = [y0,cx]
			if Voronoi.debug: print('isRightOf',self.p0,self.p1,self.direct,px,y0,cx)
			return cx>px
	def Intersect(self,Next,Voronoi_atol):
		if self.p1!=Next.p0: sys.exit('ERROR: not consecutive')
		if self.p0==(-1,-1): return None
		if Next.p1==(-2,-2): return None
		if self.direct<=0 and Next.direct>=0: return None
		if self.p0[0]==self.p1[0]==Next.p1[0]: return None
		#1: a<0 b>0  2: a<0 b=0  3: a=0 b>0  4: a=0 b=0
		ux,uy = self.p0
		vx,vy = self.p1
		wx,wy = Next.p1
		if uy>wy: ux,uy,wx,wy=wx,wy,ux,uy
		if vy>wy: vx,vy,wx,wy=wx,wy,vx,vy
		if uy>vy: ux,uy,vx,vy=vx,vy,ux,uy
		if Voronoi.debug: print('Intersect',(ux,uy),(vx,vy),(wx,wy))
		cx,ytop,cy=CalIntersect(ux,uy,vx,vy,wx,wy)
		if cx==-np.inf and ytop==np.inf and cy==np.inf: return None
		if self.direct==1  and cx+Voronoi_atol<self.base[0] \
		or Next.direct==1  and cx+Voronoi_atol<Next.base[0] \
		or self.direct==-1 and cx-Voronoi_atol>self.base[0] \
		or Next.direct==-1 and cx-Voronoi_atol>Next.base[0]:
			return None
		if self.direct==0  and cy+Voronoi_atol<self.base[1] \
		or Next.direct==0  and cy+Voronoi_atol<Next.base[1]:
			return None
		return cx,ytop

	def oldIntersect(self,Next):
		if self.p1!=Next.p0: sys.exit('ERROR: not consecutive')
		if self.p0==(-1,-1): return None
		if Next.p1==(-2,-2): return None
		if self.direct<=0 and Next.direct>=0: return None #not Intersect
		if self.p0[0]==self.p1[0]==Next.p1[0]:
			return None
			#if self.direct==Next.direct==1: return np.inf,(self.p0[1]+self.p1[1]+Next.p0[1]+Next.p1[1])/4
			#elif self.direct==Next.direct==-1: return -np.inf,(self.p0[1]+self.p1[1]+Next.p0[1]+Next.p1[1])/4
		#1: a<0 b>0  2: a<0 b=0  3: a=0 b>0  4: a=0 b=0
		ux,uy = self.p0
		vx,vy = self.p1
		wx,wy = Next.p1
		if uy>wy: ux,uy,wx,wy=wx,wy,ux,uy
		if vy>wy: vx,vy,wx,wy=wx,wy,vx,vy
		if uy>vy: ux,uy,vx,vy=vx,vy,ux,uy
		if Voronoi.debug: print('Intersect',(ux,uy),(vx,vy),(wx,wy))

		if uy == vy and ux != vx and vy != wy:
			cx = (ux + vx)/2
			if wx==cx:
				#print('1')
				cy = (uy+wy)/2 - (ux-vx)**2/8./(wy-uy)
				ytop=wy
			else:
				#solve[(x - u1)^2+(y - u2)^2=r2,(x-v1)^2+(y-u2)^2=r2, (x-w1)^2+(y-w2)^2=r2, x,y]
				#cy = (-ux*wx - vx*wx + wx**2 + wy**2 + ux*vx - uy**2)/2./(wy-uy)
				#print('2')
				cy = (wy+uy)/2. +  (ux-wx)*(vx-wx)/2./(wy-uy)
				r = np.sqrt(((ux-wx)**2+(uy-wy)**2)*((wx-vx)**2+(uy-wy)**2))/2./(wy-uy)
				ytop=cy+r
		elif wx==ux and wx!=vx:
			cy = (wy+uy)/2
			if vy==cy:
				cx = (ux+vx)/2 + (uy-wy)**2/(ux-vx)/8.
				r = ((uy-wy)**2+4*(ux-vx)**2)/8./abs(ux-vx)
				#print('4')
				ytop=cy+r
			else:
				cx = (ux+vx)/2 + (vy-wy)*(uy-vy)/(ux-vx)/2 # 2(x-(ux+vx)/2) / (vy-wy) = (uy-vy) / (ux-vx) #same angle
				r = np.sqrt(((ux-vx)**2+(uy-vy)**2)*((wy-vy)**2+(ux-vx)**2))/2./abs(ux-vx)
				#print('5')
				ytop=cy+r
		else:
			a = (ux-vx)*(vy-wy)-(uy-vy)*(vx-wx)
			#print('lewton',a)
			#if a<=0: return None #!!!!!!!!!!!!!!check
			if a==0: return None #!!!!!!!!!!!!!!check
			b1 = (ux-vx)*(ux+vx)+(uy-vy)*(uy+vy)
			b2 = (vx-wx)*(vx+wx)+(vy-wy)*(vy+wy)
			cx = (b1*(vy-wy)-b2*(uy-vy))/2./a
			cy = (b2*(ux-vx)-b1*(vx-wx))/2./a
			r = np.sqrt((cx-ux)**2+(cy-uy)**2)
			ytop=cy+r
			#print('3')
		#if self.direct==Next.direct==1 and cx<=min(self.base[0],Next.base[0]) \
		#or self.direct==Next.direct==-1 and cx>=max(self.base[0],Next.base[0]):
		#	return None
		if self.direct==1 and round(cx,Voronoi.SLVround)<round(self.base[0],Voronoi.SLVround) \
		or Next.direct==1 and round(cx,Voronoi.SLVround)<round(Next.base[0],Voronoi.SLVround) \
		or self.direct==-1 and round(cx,Voronoi.SLVround)>round(self.base[0],Voronoi.SLVround) \
		or Next.direct==-1 and round(cx,Voronoi.SLVround)>round(Next.base[0],Voronoi.SLVround):
			#print(color('Nothing11',31,1),cx,cy,self.base,Next.base)
			return None
		if self.direct==0 and round(cy,Voronoi.SLVround)<round(self.base[1],Voronoi.SLVround) \
		or Next.direct==0 and round(cy,Voronoi.SLVround)<round(Next.base[1],Voronoi.SLVround):
			#print(color('Nothing22',31,1),cx,cy,self.base,Next.base)
			return None
		return cx,ytop
	#END_OF_def oldIntersect(self,Next)
	def complete(self,xs,ys,echo=False):
		#fill self.summit
		xb,yb = self.base
		#if xs>xb and not self.direct>0: print self,xs,ys
		#if xs==xb and not self.direct==0: print self,xs,ys
		#if xs<xb and not self.direct<0: print self,xs,ys
		#The directs are good due to the checking here
		if xb<=-0.5 and xs<=-0.5 or xb>=self.RightLimit and xs>=self.RightLimit: self.summit = None
		elif yb<=-0.5 and ys<=-0.5 or yb>=self.TopLimit and ys>=self.TopLimit: self.summit = None
		#elif round(xs,Voronoi.SLVround)==round(xb,Voronoi.SLVround) and round(ys,Voronoi.SLVround)==round(yb,Voronoi.SLVround): self.summit = None #seems useless
		else:
			if not -0.5<=xb<=self.RightLimit:
				print(xb,yb,xs,ys,'Unexpected base x',xb) #I don't know why x is so good
			if yb<-0.5 and ys>self.TopLimit: #crossing the whole image
				#if echo: print(self.base,xs,ys,'->', end=' ')
				self.base = xb-(yb+0.5)*(xs-xb)/(ys-yb), -0.5
				xs,ys =  xb+(self.TopLimit-yb)*(xs-xb)/(ys-yb), self.TopLimit
				xb,yb = self.base
				#if echo: print(self.base,xs,ys)
			elif ys<-0.5 and yb>self.TopLimit: #backward crossing the whole image
				xs,ys =  xb-(yb+0.5)*(xs-xb)/(ys-yb), -0.5
				self.base = xb+(self.TopLimit-yb)*(xs-xb)/(ys-yb), self.TopLimit
				xb,yb = self.base
			elif ys>-0.5 and yb<-0.5: #upward crossing y==0 (bottom)
				if echo: print('upward')
				if self.direct==0:
					self.base = xb,-0.5
				else:
					self.base = xb-(yb+0.5)*(xs-xb)/(ys-yb), -0.5
				xb,yb = self.base
			elif yb<self.TopLimit and ys>self.TopLimit: #crossing y==max
				if echo: print(xb,yb,xs,ys,'T->')
				if self.direct==0:
					xs,ys = xb,self.TopLimit
				else:
					xs,ys = xb+(self.TopLimit-yb)*(xs-xb)/(ys-yb), self.TopLimit
				if echo: print(xb,yb,xs,ys)
			elif ys<-0.5 and yb>-0.5: #back crossing y==0
				xs,ys =  xb-(yb+0.5)*(xs-xb)/(ys-yb), -0.5
				if echo: print('back crossing')
			elif yb>self.TopLimit and ys<self.TopLimit: #back crossing y==max
				if echo: print(xb,yb,xs,ys,'T<-')
				self.base = xb+(self.TopLimit-yb)*(xs-xb)/(ys-yb), self.TopLimit
				xb,yb = self.base
				if echo: print(xb,yb,xs,ys)
			else:
				if echo: print('other1',xb,yb,xs,ys)
			#treat x y in two steps, to cover such cases that crossing two image edges
			if xs<-0.5 and xb>-0.5:
				xs,ys =  -0.5, yb-(xb+0.5)*(ys-yb)/(xs-xb)
			elif xs>self.RightLimit and xb<self.RightLimit:
				#print(xb,yb,xs,ys,'R->',end='')
				xs,ys =  self.RightLimit, yb+(self.RightLimit-xb)*(ys-yb)/(xs-xb)
				#print(xb,yb,xs,ys)
			else:
				if echo: print('other',xb,yb,xs,ys)
				if not (-0.5<=yb<=self.TopLimit and -0.5<=xb<=self.RightLimit): print('What',xb,yb,xs,ys)
			self.summit = xs,ys
	#END_OF_def complete(self,xs,ys,echo=False):

class SweepTable(dict):
	#a dict of half edges
	def __init__(self,Px,Py,Right,Top,atol=1e-8):
		self.KL = [] #key list
		self.ysweep=0
		self.Voronoi_atol=atol
		self.p = None
		self.RightLimit=Right
		self.TopLimit=Top
		self.HalfEdge=type('HalfEdge',(HalfEdgeM,),dict(RightLimit=self.RightLimit,TopLimit=self.TopLimit))

		ymin = Py>Py.min()
		self.Q = EventQueue(Px[ymin],Py[ymin],10**-Voronoi.SLVround)

		#initialize self dict
		ymin = list(np.where(Py==Py.min())[0])
		if len(ymin)==1: n=ymin.pop(0)
		else: n=ymin.pop(np.argmin(Px[ymin]))
		x0,y0 = Px[n],Py[n]
		self[-1] = self.HalfEdge((-1,-1),(x0,y0),-1,(-1,-1)) #0: left boundary of image
		#print('HalfEdge',-1,(-1,-1),(x0,y0),-1,(-1,-1))
		self.KL.append(-1)
		self.counter = itertools.count(1) #to make sure each key is different
		for n in ymin:
			x = Px[n]
			y = Py[n]
			count = next(self.counter)
			self[count] = self.HalfEdge((x0,y0),(x,y),0,((x0+x)/2.,-np.inf))
			#print('HalfEdge',count,(x0,y0),(x,y),0,((x0+x)/2.,-0.5))
			self.KL.append(count)
			x0,y0 = x,y
		self[-2] = self.HalfEdge((x0,y0),(-2,-2),1,(-2,-2)) #1: right boundary of image
		#print('HalfEdge',-2,(x0,y0),(-2,-2),1,(-2,-2))
		self.KL.append(-2)

	def __str__(self):
		s = "SweepTableStart------------------------\n"
		for k in self.KL:
			s += "\t"+str(k)+' : '+self[k].__str__()+'\n'
		s += "SweepTableEnd-------------------------"
		return s

	def locate(self,p):
		#p is a Site event
		#find the key of first edge behind which p locates
		#X(key) <= px < X(next key)
		Nlow = 0
		Nhigh = len(self.KL)-1
		#print Nlow,'--',Nhigh
		while True:
			if Nlow+1 == Nhigh:
				break
			Nmid = (Nlow+Nhigh)//2
			#print 'mid:',Nmid
			if self[self.KL[Nmid]].isRightOf(p[0],p[1],self.Voronoi_atol):
				Nhigh = Nmid
			else:
				Nlow = Nmid
		return Nlow
		'''
		#NO, seems wrong
		if (self[self.KL[Nlow+1]].p0,self[self.KL[Nlow+1]].p1)==(self[self.KL[Nlow]].p1,self[self.KL[Nlow]].p0):
			assert (self[self.KL[Nlow+1]].base,self[self.KL[Nlow+1]].direct)==(self[self.KL[Nlow]].base,-self[self.KL[Nlow]].direct)
			return Nlow+1
		else: return Nlow
		'''
	#END_OF_def locate(self,p):

	def NewEdges(self,n,x,y,otherside=1):
		#input: self[KL[n]].X <= p < self[KL[n+1]].X
		#p: y,x,...
		x0,y0 = self[self.KL[n]].p1
		base = x,(y**2-y0**2-(x-x0)**2)/2./(y-y0)
		count = next(self.counter)
		self.KL.insert(n+1,count)
		self[count] = self.HalfEdge((x0,y0),(x,y),-otherside,base,count+1) #new minus half edge
		if Voronoi.debug: print('new: ',count,self[count])
		count = next(self.counter)
		self.KL.insert(n+2,count)
		self[count] = self.HalfEdge((x,y),(x0,y0),otherside,base,count-1)  #new plus half edge
		if Voronoi.debug: print('new: ',count,self[count])

	@staticmethod
	def merge2one(HE_l,HE_r,x,y):
		#two halfedge intersect and merge into one at the Vertex event x,y
		#(x,y) is the top point of the circle which goes through the three points (l.p0,l.p1,r.p1), l.p1==r.p0
		#the center of the circle is (x,y-r)
		#(x-x0)^2 + (y-r-y0)^2 = r^2, where x0,y0 is any one of the three points
		assert HE_l.p1==HE_r.p0
		if HE_l.p0[1]==HE_r.p1[1]:
			direct = 0
		else:
			if HE_l.p0[1]>HE_r.p1[1]:
				x0 = HE_l.p0[0]
				if x==x0: direct = 1
			else:
				x0 = HE_r.p1[0]
				if x==x0: direct = -1
			if x < x0: direct = -1
			elif x > x0: direct = 1
		x0,y0 = HE_l.p0
		if abs(y-HE_l.p1[1])>abs(y-y0): x0,y0=HE_l.p1
		if abs(y-HE_r.p1[1])>abs(y-y0): x0,y0=HE_r.p1
		xb,yb = x,(y**2-y0**2-(x-x0)**2)/2./(y-y0)
		#print(HE_l)
		#print(HE_r)
		#print(x,y,x0,y0,'-->',xb,yb)
		return direct,xb,yb

	#def shrink(self): x,y,l,r = self.p n = self.KL.index(l) if r!=self.KL[n+1]: print(self.KL) print(self.KL[n],self[self.KL[n]]) print(self.KL[n+1],self[self.KL[n+1]]) print(self.KL[n+2],self[self.KL[n+2]]) assert r==self.KL[n+1] L = self.KL[n-1] R = self.KL[n+2] #print(L,l,r,R,self) #self.Q.deleteifhas(L) #self.Q.deleteifhas(r) direct,xb,yb=SweepTable.merge2one(self[l], self[r],x,y) if self[l].direct==self[r].direct==1 and xb<=min(self[l].base[0],self[r].base[0]) \ or self[l].direct==self[r].direct==-1 and xb>=max(self[l].base[0],self[r].base[0]): return [] if Voronoi.debug: print('merge:',l,r,'->',(xb,yb)) count = next(self.counter) self[count] = self.HalfEdge(self[l].p0,self[r].p1,direct,(xb,yb)) #print('HalfEdge',self[count]) self.KL[n] = count self.KL.pop(n+1) #print('comp',self[l],xb,yb) self[l].complete(xb,yb) #print('comp',self[r],xb,yb) self[r].complete(xb,yb) i1 = self[L].Intersect(self[count]) if i1 is not None: if round(i1[1],Voronoi.SLVround)<=round(self.p[1],Voronoi.SLVround): warnings.warn(f'WARNING: very small step {self.p} {i1}') if Voronoi.debug: print('->V L',i1,self[L],self[count]) self.Q.insert(i1[0],i1[1],L,count) i2 = self[count].Intersect(self[R]) if i2 is not None: if round(i2[1],Voronoi.SLVround)<=round(self.p[1],Voronoi.SLVround): warnings.warn(f'WARNING: very small step {self.p} {i2}') if Voronoi.debug: print('->V R',i2,self[count],self[R]) self.Q.insert(i2[0],i2[1],count,R) return l,r 

	#def RenewSideBoundary(self): y = self.p[1] #print(color('lt',31,1),self[self.KL[-2]],self.RightLimit,y,self.ysweep) todo = [] if y <= self.ysweep or len(self)<=2: return todo if self[self.KL[-2]].isRightOf(self.RightLimit,y): #if Voronoi.debug: print(color('Edge:',31,1),self[-2],self[self.KL[-2]]) assert self[-2].p0==self[self.KL[-2]].p1 u0,u1 = self[self.KL[-2]].p0 v0,v1 = self[self.KL[-2]].p1 assert u1!=v1 yb = (u0**2-2*u0*(self.RightLimit)+u1**2-v0**2+2*v0*(self.RightLimit)-v1**2)/2./(u1-v1) self[self.KL[-2]].complete(self.RightLimit,yb) todo.append([self.KL[-2],self[self.KL[-2]].copy()]) #maybe the VertexEvt to delete is just the p which has just been deleted from Q as the minimum. if self.p is not None and len(self.p)==4 and self.KL[-3]==self.p[2]: self.p = None print('WARNING no',self.p) #else: self.Q.deleteifhas(self.KL[-3]) self[-2] = self.pop(self.KL[-2]) self.KL.pop(-2) self[-2].p1 = (-2,-2) self[-2].direct = 1 if Voronoi.debug: print('-->',self[-2]) if not self[self.KL[1]].isRightOf(-0.5,y): #if Voronoi.debug: print(color('Edge:',31,1),self[-1],self[self.KL[1]],-0.5,y) assert self[-1].p1==self[self.KL[1]].p0 u0,u1 = self[self.KL[1]].p0 v0,v1 = self[self.KL[1]].p1 assert u1!=v1 yb = (u0**2-2*u0*(-0.5)+u1**2-v0**2+2*v0*(-0.5)-v1**2)/2./(u1-v1) self[self.KL[1]].complete(-0.5,yb) todo.append([self.KL[1],self[self.KL[1]].copy()]) if self.p is not None and len(self.p)==4 and self.KL[1]==self.p[2]: self.p = None print('WARNING no',self.p) #else: self.Q.deleteifhas(self.KL[1]) self[-1] = self.pop(self.KL[1]) self.KL.pop(1) self[-1].p0 = (-1,-1) self[-1].direct = -1 #if Voronoi.debug: print('-->:',self[-1]) self.ysweep = y return todo #END_OF_def RenewSideBoundary(self):

	#def RenewRightBoundary(self): #y = self.p[1] #print(color('lt',31,1),self[self.KL[-2]],self.RightLimit,y,self.ysweep) #todo = [] #if y <= self.ysweep or len(self)<=2: return todo #if self[self.KL[-2]].isRightOf(self.RightLimit,y): if 1: if Voronoi.debug: print(color('EdgeR:',31,1),self[-2],self[self.KL[-2]]) assert self[-2].p0==self[self.KL[-2]].p1 u0,u1 = self[self.KL[-2]].p0 v0,v1 = self[self.KL[-2]].p1 assert u1!=v1 yb = (u0**2-2*u0*(self.RightLimit)+u1**2-v0**2+2*v0*(self.RightLimit)-v1**2)/2./(u1-v1) self[self.KL[-2]].complete(self.RightLimit,yb) todo=[self.KL[-2],self[self.KL[-2]].copy()] self[-2] = self.pop(self.KL[-2]) self.KL.pop(-2) self[-2].p1 = (-2,-2) self[-2].direct = 1 if Voronoi.debug: print('-->',self[-2]) return todo #END_OF_def RenewRightBoundary(self):
	#def RenewLeftBoundary(self): y = self.p[1] #print(color('lt',31,1),self[self.KL[-2]],self.RightLimit,y,self.ysweep) #todo = [] #if y <= self.ysweep or len(self)<=2: return todo #if not self[self.KL[1]].isRightOf(-0.5,y): if 1: if Voronoi.debug: print(color('EdgeL:',31,1),self[-1],self[self.KL[1]],-0.5,y) assert self[-1].p1==self[self.KL[1]].p0 u0,u1 = self[self.KL[1]].p0 v0,v1 = self[self.KL[1]].p1 assert u1!=v1 yb = (u0**2-2*u0*(-0.5)+u1**2-v0**2+2*v0*(-0.5)-v1**2)/2./(u1-v1) self[self.KL[1]].complete(-0.5,yb) todo=[self.KL[1],self[self.KL[1]].copy()] #if self.p is not None and len(self.p)==4 and self.KL[1]==self.p[2]: self.p = None print('WARNING no',self.p) #else: self.Q.deleteifhas(self.KL[1]) self[-1] = self.pop(self.KL[1]) self.KL.pop(1) self[-1].p0 = (-1,-1) self[-1].direct = -1 #if Voronoi.debug: print('-->:',self[-1]) return todo #END_OF_def RenewLeftBoundary(self):

class EventQueue(object):
	def __init__(self,Px,Py,Voronoi_atol=1e-8):
		if Px.shape==Py.shape==(Px.size,):
			self.S = np.zeros((Px.size,2))
			self.Stotal = Px.size
			self.S[:,0] = Py
			self.S[:,1] = Px
			self.SMin = 0
			self.counter = itertools.count(1)
			self.NumEvt = 0 #?useful?
			self.VerEvt = []
			self.Voronoi_atol=Voronoi_atol
		else:
			sys.exit('ERROR: shape error')
	def __str__(self):
		return str(self.VerEvt)
		#return self.__repr__()+'\n'+str(self.VerEvt)
	def delMin(self):
		while self.NumEvt>0 and self.VerEvt[0][2] == 0:
			heapq.heappop(self.VerEvt)
			self.NumEvt-=1
		#The interesting property of a heap is that its smallest element is always the root, heap[0].
		if self.SMin<self.S.shape[0] and (self.NumEvt==0 or self.S[self.SMin].tolist() <= self.VerEvt[0][:2]):
			self.SMin += 1
			return self.S[self.SMin-1][[1,0]]	#Site event
		elif self.NumEvt>0 and (self.SMin>=self.S.shape[0] or self.S[self.SMin].tolist() > self.VerEvt[0][:2]):
			e = heapq.heappop(self.VerEvt)
			self.NumEvt-=1
			if self.NumEvt>0 and e[0]-self.Voronoi_atol>=self.VerEvt[0][0]: warnings.warn(f'WARNING: next step is very near {e} {self.VerEvt[0]}')
			#self.pop(n) #!!! seems useless
			ne=len(e)
			if ne==5:
				y,x,count,n,nr=e
				return [x,y,n,nr]	#Vertex event
			elif ne==6: #turns out to be useless
				y,x,count,n,nr,nl=e
				return [x,y,n,nr,nl]
			else: sys.exit('ERR')
		else:
			return None
		#END_OF_def delMin(self):
	def insert(self,x,y,n,nr): #simplified again
		count=next(self.counter)
		heapq.heappush(self.VerEvt,[y,x,-count,n,nr]) #???-count
		self.NumEvt+=1
	#def insert(self,x,y,n,nr,nl=None): #four or five elements #turns out this is useless
	#	count=next(self.counter)
	#	if nl is None: heapq.heappush(self.VerEvt,[y,x,count,n,nr])
	#	else: heapq.heappush(self.VerEvt,[y,x,count,n,nr,nl])
	#def insert(self,x,y,l,r,count=None): #four elements #add an option of five
	#	#assert l not in self #!!! seems useless
	#	if count is None: count=next(self.counter)
	#	evt = [y,x,count,l,r]
	#	#self[l] = evt #!!! seems useless
	#	heapq.heappush(self.VerEvt,evt)
	#def deleteifhas(self,l):
		#pass #!!! seems useless
		#evt = self.pop(l,None)
		#if evt is not None: evt[2] = 0
#END_OF_class EventQueue(dict):

class Voronoi(object):
	SLVround = 6
	#if SLVround is low, two points far from each other can be taken as neighbors
	debug = False
	#debug = True
	def __init__(self,image=None,events=None,**kwargs):
		"""
		image: a 2D array in which >0 means non-empty.
		events: coordinates of points, 2D array with shape (# of points, 2).
		if calArea=True, calculate area of each cell.
		if calCentroid=True, calculate the centroid of each cell.
		Resolution: int, number of decimals to round on the input event coordinates
			default=3. But no rounding on input when "calCentroid".
		"""
		#StartTime=time.time()
		if kwargs.get('Silent',True): warnings.simplefilter('ignore')
		self.FileName = kwargs.pop('FileName','tmp')
		self.ToCalArea = kwargs.pop('calArea',False)
		self.ToCalPVD = False #kwargs.pop('calpvd',False)
		self.ToCalDVD = kwargs.pop('caldvd',False)
		self.ToCalAdj = kwargs.pop('caladj',False)
		self.ToCalDst = False #kwargs.get('caldst',False)
		self.ToCalTri = kwargs.get('calTriangle',False)
		self.ToCalCtd = kwargs.get('calCentroid',False)
		self.ToCum2D = kwargs.get('cum2d',False)
		self.ToCalDel = kwargs.pop('Delaunay',False)
		self.SmoothNumber = kwargs.get('SmoothNumber',0)
		self.SmoothFactor = kwargs.get('SmoothFactor',0.2)
		self.ToRemoveEdgePoint = kwargs.get('RemoveEdgePoint',False)
		self.ToMarkEdgePoint = kwargs.get('MarkEdgePoint',False)
		self.ToGetPoly = kwargs.get('GetPoly',False)
		self.PolyFile = kwargs.pop('PolyFile',None)
		self.Hdr = kwargs.get('Hdr',None)
		self.OffSetX=0
		self.OffSetY=0
		self.scale=1
		self.Voronoi_atol=10**-Voronoi.SLVround
		if image is not None: #integer position
			self.Mode='image'
			assert events is None
			if kwargs.pop('CleanFilledBorder',False):
				n0=np.where(np.sum(image>0,axis=0)<image.shape[0])[0]
				n1=np.where(np.sum(image>0,axis=1)<image.shape[1])[0]
				if n1[0]-1>0: image[:n1[0]-1,:]=0
				if n1[-1]+2<image.shape[0]-1: image[n1[-1]+2:,:]=0
				if n0[0]-1>0: image[:,:n0[0]-1]=0
				if n0[-1]+2<image.shape[1]-1: image[:,n0[-1]+2:]=0
			self.ImgSizeX,self.ImgSizeY = image.shape
			self.RightLimit=self.ImgSizeX-0.5
			self.TopLimit=self.ImgSizeY-0.5
			self.pmap = np.zeros(image.shape,dtype='int32') #pvd
			self.Py,self.Px = np.where(image.T>0) #coordinates (y,x)
			self.pmap[self.Px,self.Py] = np.arange(1,self.Px.size+1)
		elif events is not None: #float position
			self.Mode='event'
			assert events.shape[1] >= 2
			#events=events[:,:2] #2022.05.07
			if self.ToCalPVD: sys.exit('--calpvd not supported in case of events (as opposed to image) input')
			if self.ToCalDst: sys.exit('--caldst not supported in case of events (as opposed to image) input')

			border=kwargs.get('border',None)
			if border is None:
				if kwargs.get('autoscale',True):
					self.scale = 2*np.sqrt(events.shape[0])/(np.max(events[:,:2])-np.min(events[:,:2]))
					if self.scale<2: self.scale=1
					else:
						self.scale=np.round(self.scale,1)
						#print(f'Scale the coordinates by: {self.scale}')
						events[:,:2] *= self.scale
				xlow = np.floor(np.min(events[:,0]))-1
				ylow = np.floor(np.min(events[:,1]))-1
				xhigh = np.ceil(np.max(events[:,0]))+1
				yhigh = np.ceil(np.max(events[:,1]))+1
				self.OffSetX = -0.5-xlow #shift the lower border to -0.5
				events[:,0] += self.OffSetX
				xhigh += self.OffSetX
				self.OffSetY = -0.5-ylow #shift the lower border to -0.5
				events[:,1] += self.OffSetY
				yhigh += self.OffSetY
				self.ImgSizeX = int(xhigh+0.5)
				self.ImgSizeY = int(yhigh+0.5)
				self.RightLimit=self.ImgSizeX-0.5
				self.TopLimit=self.ImgSizeY-0.5
				#np.savetxt('try',events)
				#https://docs.python.org/3/tutorial/floatingpoint.html
				#Some float numbers do not have exact representations in binary floating point, e.g., float.hex(0.1): '0x1.999999999999ap-4' 0.1!=1-0.9. Luckily int-0.5 seems OK, float.hex(0.5): '0x1.0000000000000p-1'.
			else:
				xlow = border['xlow']
				ylow = border['ylow']
				xhigh = border['xhigh']
				yhigh = border['yhigh']
				self.OffSetX = -0.5-xlow #shift the lower border to -0.5
				events[:,0] += self.OffSetX
				xlow += self.OffSetX
				xhigh += self.OffSetX
				self.OffSetY = -0.5-ylow #shift the lower border to -0.5
				events[:,1] += self.OffSetY
				ylow += self.OffSetY
				yhigh += self.OffSetY
				self.ImgSizeX = int(xhigh)+1
				self.ImgSizeY = int(yhigh)+1
				self.RightLimit = xhigh
				self.TopLimit = yhigh
				#print(f'Set image size: {xlow:g}~{xhigh:g} {ylow:g}~{yhigh:g}')
			print(color('Voronoi Tessellation: '+self.FileName,34,1),end=' ')
			print(f"OffSet: {self.OffSetX} {self.OffSetY} Rescale: {self.scale} ImgSize: {self.ImgSizeX} {self.ImgSizeY}")
			if np.min(events[:,0])<-0.5 or np.min(events[:,1])<-0.5 or np.max(events[:,0])>self.RightLimit or np.max(events[:,1])>self.TopLimit:
				print(color(f"ERROR: points out of -0.5~{self.RightLimit:g}, -0.5~{self.TopLimit:g}",31,1))
				exit()

			self.MakeIntImage=kwargs.pop('MakeIntImage',-1) #-1 means do nothing
			if self.MakeIntImage>=0: #can be 0 (python) or 1 (ds9)
				if self.ImgSizeX>4096 or self.ImgSizeY>4096:
					warnings.warn("The image size exceeds 4096. Do you really want to make such a large image?")
					#readyn=input()
					#if readyn!='y': exit()
				if self.ImgSizeX<3 or self.ImgSizeY<3:
					warnings.warn("The image size < 3. Do you really want to make such a small image?")
				pmap = np.zeros((int(self.ImgSizeX),int(self.ImgSizeY)),dtype='float64')
				if Voronoi.debug: print("init pmap",pmap.shape)
				with open(self.FileName+'_points.reg','w') as freg:
					for x,y in events[:,:2]:
						if self.MakeIntImage==0: 
							pmap[int(x+0.5),int(y+0.5)] += 1
							print(f"image;circle({y+1:.2f},{x+1:.2f},0.5) # color=red", file=freg)
						elif self.MakeIntImage==1:
							pmap[int(y+0.5),int(x+0.5)] += 1
							print(f"image;circle({x+1:.2f},{y+1:.2f},0.5) # color=red", file=freg)
				if self.Hdr is None: self.Hdr=fits.Header()
				self.Hdr.update({'OffSetX':self.OffSetX})
				self.Hdr.update({'OffSetY':self.OffSetY})
				self.Hdr.update({'Scale':self.scale})
				fits.writeto(self.FileName+'_image.fits',pmap,self.Hdr,overwrite=True)
				print('>> '+self.FileName+'_image.fits')
				print('>> '+self.FileName+'_points.reg')
			Resolution=kwargs.pop('Resolution',3)
			if not self.ToCalCtd:
				events[:,:2]=np.round(events[:,:2],Resolution)
			tmp,ind,cts=np.unique(events[:,1::-1],return_index=True,return_counts=True,axis=0) #switch xy
			if tmp.shape[0]!=events.shape[0]:
				warnings.warn(color("found %d duplicated points" % (events.shape[0]-len(tmp)),31,1))
			self.Py,self.Px = tmp[:,0],tmp[:,1] #switch back
			#if events.shape[1] == 3: value = tmp[:,2] #to be continued
			self.pmap={(self.Px[n],self.Py[n]):n+1 for n in np.arange(self.Px.size)}
			self.pind={(self.Px[n],self.Py[n]):[ind[n],cts[n]] for n in np.arange(self.Px.size)}
			if events.shape[1]>2: self.pID=events[ind,2]
			else: self.pID=None
		else: sys.exit('ERROR: image mode or events mode?')

		self.Edges = {}
		self.Amap = None
		self.Wmap = None
		self.Smap = None
		self.Vmap = None
		self.filledmap=None
		self.ctdvector=None
		self.PPdis = None
		self.Adj = None
		self.EdgePoint = None
		if self.ToCalDst: #image mode, not events mode
			self.ToCalArea = True # need self.Amap
			self.ToCalPVD = True # need PVD to fill all the cells in the image
		if self.ToCum2D or self.ToCalTri or self.ToCalCtd or self.SmoothNumber>0:
			self.ToCalArea = True # need self.Amap
		if self.ToCalDVD:
			self.ToCalArea = True
			self.ToGetPoly = True

		self.Construct(**kwargs)
		if self.ToCalArea: self.CalArea(**kwargs)
		if self.SmoothNumber>0: self.SmoothDensity(**kwargs)
		if self.ToCalPVD: self.CalPVD(**kwargs)
		if self.ToCalDVD: self.CalDVD(**kwargs)
		if self.ToCalAdj: self.CalAdj()
	#END_OF_init

	def saveresults(self):
		if self.Mode=='image' or self.MakeIntImage>=0:
			d = np.array([[e.base[1],e.base[0],e.summit[1],e.summit[0],e.p0[1],e.p0[0],e.p1[1],e.p1[0]] for e in self.Edges.values() if e.summit is not None])+1 #-np.array([self.OffSetY,self.OffSetX,self.OffSetY,self.OffSetX,self.OffSetY,self.OffSetX,self.OffSetY,self.OffSetX])
			open(self.FileName+'.reg','w').write('image\n')
			np.savetxt(self.FileName+'.reg',d,fmt="line(%.3f,%.3f,%.3f,%.3f) # tag={%g,%g,%g,%g}",header='image',comments='')
			print('>> '+self.FileName+'.reg')
			if self.ToCalDel and False:
				d = np.array([[e.p0[1],e.p0[0],e.p1[1],e.p1[0],e.base[1],e.base[0],e.summit[1],e.summit[0]] for e in self.Edges.values() if e.summit is not None])+1-np.array([self.OffSetY,self.OffSetX,self.OffSetY,self.OffSetX,self.OffSetY,self.OffSetX,self.OffSetY,self.OffSetX])
				open(self.FileName+'_Delaunay.reg','w').write('image\n')
				np.savetxt(self.FileName+'_Delaunay.reg',d,fmt="line(%g,%g,%g,%g) # tag={%.3f,%.3f,%.3f,%.3f}")
				print('>> '+self.FileName+'_Delaunay.reg')
		if self.Mode=='event' or self.Mode=='image':
			with open(self.FileName+'_VT.dat','w') as fout:
				print('#'+json.dumps({'xlow':(-0.5-self.OffSetX)/self.scale,'ylow':(-0.5-self.OffSetY)/self.scale,'xhigh':(self.RightLimit-self.OffSetX)/self.scale,'yhigh':(self.TopLimit-self.OffSetY)/self.scale,'scale':self.scale}),file=fout)
				#for debug: d=np.array([[e.base[0],e.base[1],e.summit[0],e.summit[1],e.p0[0],e.p0[1],e.p1[0],e.p1[1]] for e in self.Edges.values() if e.summit is not None])
				d=np.array([[e.base[0],e.base[1],e.summit[0],e.summit[1],e.p0[0],e.p0[1],e.p1[0],e.p1[1]] for e in self.Edges.values() if e.summit is not None])-[self.OffSetX,self.OffSetY,self.OffSetX,self.OffSetY,self.OffSetX,self.OffSetY,self.OffSetX,self.OffSetY]
				np.savetxt(fout,np.hstack((np.arange(1,len(d)+1).reshape(len(d),1),d/self.scale)),fmt="%d	%.9g %.9g %.9g %.9g	%.9g %.9g %.9g %.9g")
				print('>>',fout.name)
		if self.ToCalDst: #image mode, not events mode
			assert np.all(self.Amap[self.Px,self.Py]>0) #self.CalArea done
			assert np.all(self.pmap>0) #self.CalPVD done
			dmap = np.zeros(image.shape,dtype='float64')
			dmap[self.Px,self.Py] = image[self.Px,self.Py]/self.Amap[self.Px,self.Py]
			dmap[:,:] = dmap[self.Px,self.Py][self.pmap[self.pmap>0]-1].reshape(self.ImgSizeX,self.ImgSizeY)
			#vmap is not uniquely defined, thus by now the uncertainty is transfered into dmap
			print('>> '+self.FileName+'_dst.fits')
			fits.writeto(self.FileName+'_dst.fits',dmap,self.Hdr,overwrite=True)
		if self.ToCalPVD:
			print('>> '+self.FileName+'_pvd.fits')
			fits.writeto(self.FileName+'_pvd.fits',self.pmap,self.Hdr,overwrite=True)
		if self.ToCalArea and not self.ToCalCtd:
			if type(self.Amap) is dict:
				with open(self.FileName+'_area.dat','w') as fout:
					for (x,y) in self.Amap:
						print(f'{self.pind[x,y][0]:d}	{(x-self.OffSetX)/self.scale:.9g} {(y-self.OffSetY)/self.scale:.9g} {self.Amap[x,y]:.9g} {self.pind[x,y][1]:d}',file=fout)
					print('>>',fout.name)
			else:
				print('>> '+self.FileName+'_area.fits')
				fits.writeto(self.FileName+'_area.fits',self.Amap,self.Hdr,overwrite=True)
		if self.SmoothNumber>0:
			if type(self.Smap) is dict:
				with open(self.FileName+'_SmoFactor.dat','w') as fout:
					for (x,y) in self.Smap:
						print(f'{self.pind[x,y][0]:d}	{(x-self.OffSetX)/self.scale:.9g} {(y-self.OffSetY)/self.scale:.9g} {self.Smap[x,y]:.9g} {self.pind[x,y][1]:d}',file=fout)
					print('>>',fout.name)
			else:
				print('>> '+self.FileName+'_SmoFactor.fits')
				fits.writeto(self.FileName+'_SmoFactor.fits',self.Smap,self.Hdr,overwrite=True)
		if self.ToCalCtd:
			with open(self.FileName+'_ctd.dat','w') as fout:
				for n in range(len(self.Px)):
					x,y=self.Px[n],self.Py[n]
					if type(self.Amap) is dict:
						print(f'{self.pind[x,y][0]:d}	{(x-self.OffSetX)/self.scale:.9g} {(y-self.OffSetY)/self.scale:.9g} {(self.Wmap[x,y][0]-self.OffSetX)/self.scale:.9g} {(self.Wmap[x,y][1]-self.OffSetY)/self.scale:.9g} {self.Amap[x,y]:.9g} {self.pind[x,y][1]:d}',file=fout)
					else:
						print(f'{self.pmap[x,y]:d}	{(x-self.OffSetX)/self.scale:.9g} {(y-self.OffSetY)/self.scale:.9g} {(self.Wmap[x,y][0]-self.OffSetX)/self.scale:.9g} {(self.Wmap[x,y][1]-self.OffSetY)/self.scale:.9g} {self.Amap[x,y]:.9g}',file=fout)
				print('>>',fout.name)
			if self.Mode=='image' or self.MakeIntImage>=0:
				#Ps=np.array(list(self.Amap.keys()))
				#np.savetxt(self.FileName+'_ctd.reg', np.vstack((Ps[:,1]+1,Ps[:,0]+1,np.sqrt(np.sum(self.ctdvector**2,axis=1)),np.arctan2(self.ctdvector[:,0],self.ctdvector[:,1])*180/np.pi)).T,fmt='# vector(%f,%f,%f,%f) vector=1 width=1')
				np.savetxt(self.FileName+'_ctd.reg', np.vstack((self.Py+1,self.Px+1,np.sqrt(np.sum(self.ctdvector**2,axis=1)),np.arctan2(self.ctdvector[:,0],self.ctdvector[:,1])*180/np.pi)).T,fmt='# vector(%f,%f,%f,%f) vector=1 width=1')
				print('>> '+self.FileName+'_ctd.reg')
		if self.ToCum2D:
				print('>> '+self.FileName+'_cum2d.dat')
				np.savetxt(self.FileName+'_cum2d.dat',np.hstack((np.array(list(self.CumWmap.values())).reshape(len(self.CumWmap),1),np.array(list(self.CumWmap.keys()))-[self.OffSetX,self.OffSetY])),fmt='%f	%f	%f')
	#END_OF_def saveresults(self):

	def getarealist(self):
		if type(self.Amap) is dict:
			return [self.Amap[x,y] for (x,y) in self.Amap],[self.pind[x,y][0] for (x,y) in self.Amap],[self.pind[x,y][1] for (x,y) in self.Amap]
		else: print('later')

	def Construct(self,**kwargs):
		if Voronoi.debug: StartTime=time.time()
		print(color('Sweepline VT Construction',32,0))
		T = SweepTable(self.Px,self.Py,self.RightLimit,self.TopLimit,atol=self.Voronoi_atol)
		Q = T.Q
		pbar = tqdm(total=Q.Stotal-1)
		T.p = Q.delMin()
		#assert len(T.Q.VerEvt)==T.Q.NumEvt #remove to speed up
		ndone=1
		while T.p is not None:
			#for k,e in T.RenewSideBoundary(): #???????
			#	self.Edges[k] = e
			if T.p is None: #p is deleted during renewing side boundary edges.
				print('lksjlkfjdslkjfdslkjfdslk')
				pass
			elif len(T.p)==2: #Site Event
				if Voronoi.debug: print('SiteEvent:',T.p)
				if Voronoi.debug: print('TS1',T)
				n = T.locate(T.p)
				L = T.KL[n]
				R = T.KL[n+1]
				T.NewEdges(n,T.p[0],T.p[1])
				if Voronoi.debug: print('TS2',T)
				l = T.KL[n+1]
				r = T.KL[n+2]
				assert R == T.KL[n+3]
				#from left to right L l r R
				#Q.deleteifhas(L)
				self.afterinsertbetween(T,L,R,l,r)
				'''
				BoundaryUpdated=False
				i1 = T[L].Intersect(T[l])
				if i1 is None or i1[0]<-0.5:
					if T.p[1]>T.ysweep and L==T.KL[1] and not T[L].isRightOf(-0.5,T.p[1]):
						self.RenewLeftBoundary(T)
						BoundaryUpdated=True
					#else: print(L,T[L])
				else:
					if T[L].direct!=0 or T[L].p0[0]+T[L].p1[0]!=2*T[l].p1[0]:
						if round(i1[1],Voronoi.SLVround)<=round(T.p[1],Voronoi.SLVround): warnings.warn(f'WARNING: very small step {T.p} {i1}')
					if Voronoi.debug: print('->V L',i1,T[L],T[l])
					Q.insert(i1[0],i1[1],L,l)
				i2 = T[r].Intersect(T[R])
				if i2 is None or i2[0]>self.RightLimit:
					#print( T.p[1],T.ysweep,R,T.KL[-2],T[R].isRightOf(self.RightLimit,T.p[1]))
					if T.p[1]>T.ysweep and R==T.KL[-2] and T[R].isRightOf(self.RightLimit,T.p[1]):
						self.RenewRightBoundary(T)
						BoundaryUpdated=True
					#else: print(R,T[R])
				else:
					if round(i2[1],Voronoi.SLVround)<=round(T.p[1],Voronoi.SLVround): warnings.warn(f'WARNING: very small step {T.p} {i2}')
					if Voronoi.debug: print('->V R',i2,T[r],T[R])
					Q.insert(i2[0],i2[1],r,R)
				if BoundaryUpdated: self.ysweep = T.p[1]
				'''

				'''
				When T[L].direct==0, the SiteEvent leads to a VertexEvent which have the same y and thus needs to be processed immediately.
				The following 5-elements method is Wrong as it adds an additional pair of SiteEvent edges.
				if T[L].direct==0: assert n>=1 T.NewEdges(n-1,T.p[0],T.p[1],otherside=-1) if Voronoi.debug: print('TS2',T) l2 = T.KL[n] r2 = T.KL[n+1] i3 = T[r2].Intersect(T[L]) if i3 is not None: if Voronoi.debug: print('->V',i3,T[r2],T[L]) if Voronoi.debug: print('->V',i1,T[L],T[l]) print(i3,i1) assert round(i3[0],Voronoi.SLVround)==round(i1[0],Voronoi.SLVround) \ and round(i3[1],Voronoi.SLVround)==round(i1[1],Voronoi.SLVround) Q.insert(i3[0],i3[1],L,l,r2)
				'''

			elif len(T.p)>=4: #Vertex Event
				if T.p[2] not in T or T.p[3] not in T:
					pass
					#print('Useless VertexEvent:',T.p)
				else:
					for n in T.p[2:]: assert n in T
					if Voronoi.debug: print('VertexEvent:',T.p[0],T.p[1])
					if Voronoi.debug: print('TV1',T.p,T)
					for n in self.shrink(T):
						#print('pop',T[n])
						self.Edges[n] = T.pop(n)
					if Voronoi.debug: print('TV2',T)
			else: sys.exit('error2')
			#if Voronoi.debug: print('T',T)
			#if Voronoi.debug: print(color('Q',31,1),Q)
			T.p = Q.delMin()
			#assert len(T.Q.VerEvt)==T.Q.NumEvt #remove to speed up
			if Q.SMin>ndone:
				if (Q.SMin-pbar.n-1)*100>Q.Stotal: pbar.update(Q.SMin-ndone)
				ndone=Q.SMin
		if pbar.total>pbar.n: pbar.update(pbar.total-pbar.n)
		pbar.close()
		#END_OF_while p is not None:
		T.pop(-1)
		T.pop(-2)
		self.finish(T)
		self.clean()
	#END_OF_Construct()

	def finish(self,T):
		#???
		Voronoi_atol=10**-Voronoi.SLVround
		#for e in self.Edges.values(): #Useless? It's done again later. Try remove it.
		#	if e.summit is not None and abs(e.summit[0]-e.base[0])<Voronoi_atol and abs(e.summit[1]-e.base[1])<Voronoi_atol: e.summit = None
		for k in self.Edges.keys():
			if self.Edges[k].summit is not None and (self.Edges[k].base[0]<-0.5 or self.Edges[k].base[0]>self.ImgSizeX-0.5 or self.Edges[k].base[1]<-0.5 or self.Edges[k].base[1]>self.ImgSizeY-0.5 or self.Edges[k].summit[0]<-0.5 or self.Edges[k].summit[0]>self.ImgSizeX-0.5 or self.Edges[k].summit[1]<-0.5 or self.Edges[k].summit[1]>self.ImgSizeY-0.5): print("ERROR",self.Edges[k])
		#???

		for l in list(T.keys()):
			e = T[l]
			#e.complete((e.p0[0]**2+e.p0[1]**2-e.p1[0]**2-e.p1[1]**2+2*(e.p1[1]-e.p0[1])*(self.TopLimit))/2/(e.p0[0]-e.p1[0]), self.TopLimit)
			if e.direct==0:
				#print(0,e.base,e.base[0],self.TopLimit)
				e.complete(e.base[0],self.TopLimit)
			elif e.direct==1:
				#print(1,e.base,self.RightLimit, (e.p0[0]**2+e.p0[1]**2-e.p1[0]**2-e.p1[1]**2+2*(e.p1[0]-e.p0[0])*(self.RightLimit))/2/(e.p0[1]-e.p1[1]))
				e.complete(self.RightLimit, (e.p0[0]**2+e.p0[1]**2-e.p1[0]**2-e.p1[1]**2+2*(e.p1[0]-e.p0[0])*(self.RightLimit))/2/(e.p0[1]-e.p1[1]))
			elif e.direct==-1:
				#print(-1,e.base,-0.5, (e.p0[0]**2 + e.p0[0] + e.p0[1]**2 - e.p1[0]**2 - e.p1[1]**2 - e.p1[0])/2/(e.p0[1]-e.p1[1]))
				e.complete(-0.5, (e.p0[0]**2 + e.p0[0] + e.p0[1]**2 - e.p1[0]**2 - e.p1[1]**2 - e.p1[0])/2/(e.p0[1]-e.p1[1])) #although mostly they cross self.ImgSizeY, there is a low possibility not. So I input the crossing point at Left or Right edge of image, let complete to recalculate the summit
			self.Edges[l] = T.pop(l)

		for e in self.Edges.values():
			if e.summit is not None:
				if abs(e.summit[0]-e.base[0])<Voronoi_atol and abs(e.summit[1]-e.base[1])<Voronoi_atol:
					e.summit = None
					continue
				if e.twin is not None:
					e2 = self.Edges[e.twin]
					if e2.summit is not None:
						assert abs(e.base[0]-e2.base[0])<Voronoi_atol \
						and abs(e.base[1]-e2.base[1])<Voronoi_atol
						assert self.Edges[e2.twin]==e
						e.base = e2.summit
						e2.summit = None
						assert e.p0==e2.p1 and e.p1==e2.p0
	def clean(self):
		#Delete nouse edges
		for k in list(self.Edges.keys()):
			if self.Edges[k].summit is None: self.Edges.pop(k)
			elif self.Edges[k].base[0]<-0.5 or self.Edges[k].base[0]>self.ImgSizeX-0.5 or self.Edges[k].base[1]<-0.5 or self.Edges[k].base[1]>self.ImgSizeY-0.5 or self.Edges[k].summit[0]<-0.5 or self.Edges[k].summit[0]>self.ImgSizeX-0.5 or self.Edges[k].summit[1]<-0.5 or self.Edges[k].summit[1]>self.ImgSizeY-0.5: print("ERROR",self.Edges[k])

		for k in list(self.Edges.keys()):
			if (self.Edges[k].base[0],self.Edges[k].summit[0]) in ((-0.5,self.RightLimit),(self.RightLimit,-0.5)) \
			or (self.Edges[k].base[1],self.Edges[k].summit[1]) in ((-0.5,self.TopLimit),(self.TopLimit,-0.5)):
				print('NOTE',self.Edges[k].p0,self.Edges[k].p1,self.Edges[k].base,self.Edges[k].summit)

	def RenewRightBoundary(self,T):
		if Voronoi.debug: print(color('EdgeR:',31,1),T[-2],T[T.KL[-2]])
		assert T[-2].p0==T[T.KL[-2]].p1
		u0,u1 = T[T.KL[-2]].p0
		v0,v1 = T[T.KL[-2]].p1
		assert u1!=v1
		yb = (u0**2-2*u0*(T.RightLimit)+u1**2-v0**2+2*v0*(T.RightLimit)-v1**2)/2./(u1-v1)
		T[T.KL[-2]].complete(T.RightLimit,yb)
		self.Edges[T.KL[-2]] = copy(T[T.KL[-2]])
		T[-2] = T.pop(T.KL[-2])
		T.KL.pop(-2)
		T[-2].p1 = (-2,-2)
		T[-2].direct = 1
		if Voronoi.debug: print('-->',T[-2])
	#END_OF_def RenewRightBoundary(self,T):
	def RenewLeftBoundary(self,T):
		if Voronoi.debug: print(color('EdgeL:',31,1),T[-1],T[T.KL[1]],-0.5)
		assert T[-1].p1==T[T.KL[1]].p0
		u0,u1 = T[T.KL[1]].p0
		v0,v1 = T[T.KL[1]].p1
		assert u1!=v1
		yb = (u0**2-2*u0*(-0.5)+u1**2-v0**2+2*v0*(-0.5)-v1**2)/2./(u1-v1)
		T[T.KL[1]].complete(-0.5,yb)
		self.Edges[T.KL[1]] = copy(T[T.KL[1]])
		T[-1] = T.pop(T.KL[1])
		T.KL.pop(1)
		T[-1].p0 = (-1,-1)
		T[-1].direct = -1
		if Voronoi.debug: print('-->:',T[-1])
	#END_OF_def RenewLeftBoundary(self,T):

	def afterinsertbetween(self,T,L,R,l,r=None):
		BoundaryUpdated=False
		if r is None: r=l
		i1 = T[L].Intersect(T[l],self.Voronoi_atol)
		#print(L,l,'->',i1)
		if i1 is None or i1[0]<-0.5:
			if T.p[1]>T.ysweep and L==T.KL[1] and not T[L].isRightOf(-0.5,T.p[1],self.Voronoi_atol):
				self.RenewLeftBoundary(T)
				BoundaryUpdated=True
			#else: print(L,T[L])
		else:
			if T[L].direct!=0 or T[L].p0[0]+T[L].p1[0]!=2*T[l].p1[0]:
				if i1[1]+self.Voronoi_atol<=T.p[1]: warnings.warn(f'WARNING: very small step {T.p} {i1}')
			if Voronoi.debug: print('->V L',i1,T[L],T[l])
			T.Q.insert(i1[0],i1[1],L,l)
		i2 = T[r].Intersect(T[R],self.Voronoi_atol)
		#print(r,R,'->',i2)
		if i2 is None or i2[0]>self.RightLimit:
			if T.p[1]>T.ysweep and R==T.KL[-2] and T[R].isRightOf(self.RightLimit,T.p[1],self.Voronoi_atol):
				self.RenewRightBoundary(T)
				BoundaryUpdated=True
			#else: print(R,T[R])
		else:
			if i2[1]+self.Voronoi_atol<=T.p[1]: warnings.warn(f'WARNING: very small step {T.p} {i2}')
			if Voronoi.debug: print('->V R',i2,T[r],T[R])
			T.Q.insert(i2[0],i2[1],r,R)
		if BoundaryUpdated: self.ysweep = T.p[1]
	#END_OF_def afterinsertbetween(self,T,L,R,l,r=None):

	def shrink(self,T):
		x,y,l,r = T.p
		#print('shrink4start',x,y,l,T[l],r,T[r])
		n = T.KL.index(l)
		if r!=T.KL[n+1]:
			print('shrink4starterr',x,y,l,T[l],r,T[r])
			print(T.KL)
			print(T.KL[n],T[T.KL[n]])
			print(T.KL[n+1],T[T.KL[n+1]])
			print(T.KL[n+2],T[T.KL[n+2]])
		assert r==T.KL[n+1]
		L = T.KL[n-1]
		R = T.KL[n+2]
		#print(L,l,r,R,T)
		#T.Q.deleteifhas(L)
		#T.Q.deleteifhas(r)
		direct,xb,yb=SweepTable.merge2one(T[l], T[r],x,y)
		if T[l].direct==1  and xb+self.Voronoi_atol<T[l].base[0] \
		or T[r].direct==1  and xb+self.Voronoi_atol<T[r].base[0] \
		or T[l].direct==-1 and xb-self.Voronoi_atol>T[l].base[0] \
		or T[r].direct==-1 and xb-self.Voronoi_atol>T[r].base[0]:
			return []
		if T[l].direct==0  and yb+self.Voronoi_atol<T[l].base[1] \
		or T[r].direct==0  and yb+self.Voronoi_atol<T[r].base[1]:
			return []
		if Voronoi.debug: print('merge:',l,r,'->',(xb,yb))
		count = next(T.counter)
		T[count] = T.HalfEdge(T[l].p0,T[r].p1,direct,(xb,yb))
		#print('HalfEdge',T[count])
		T.KL[n] = count
		T.KL.pop(n+1)
		#print('comp',T[l],xb,yb)
		T[l].complete(xb,yb)
		#print('comp',T[r],xb,yb)
		T[r].complete(xb,yb)
		self.afterinsertbetween(T,L,R,count)
		'''
		i1 = T[L].Intersect(T[count])
		if i1 is not None:
			if round(i1[1],Voronoi.SLVround)<=round(T.p[1],Voronoi.SLVround): warnings.warn(f'WARNING: very small step {T.p} {i1}')
			if Voronoi.debug: print('->V L',i1,T[L],T[count])
			T.Q.insert(i1[0],i1[1],L,count)
		i2 = T[count].Intersect(T[R])
		if i2 is not None:
			if round(i2[1],Voronoi.SLVround)<=round(T.p[1],Voronoi.SLVround): warnings.warn(f'WARNING: very small step {T.p} {i2}')
			if Voronoi.debug: print('->V R',i2,T[count],T[R])
			T.Q.insert(i2[0],i2[1],count,R)
		'''
		return l,r
	#END_OF_def shrink(self):

	def drawLine(self,e):
		if np.isinf(e.summit[0]): print(e)
		if e.direct == 0:
			if e.base[0]!=e.summit[0]: sys.exit('ERROR')
			if e.base[1]<e.summit[1]:
				y1 = e.base[1]
				y2 = e.summit[1]
			else:
				y2 = e.base[1]
				y1 = e.summit[1]
			x = round(e.base[0],Voronoi.SLVround)
			if x<0: x=0
			if x == int(x):
				y1 = round(y1,Voronoi.SLVround)
				if y1<0: y1=0
				elif y1==int(y1): y1 = int(y1)
				else: y1 = int(y1)+1
				y2 = int(y2)+1
				self.pmap[e.base[0],y1:y2] = -1
			return
		elif e.direct < 0:
			if e.p0[1]>e.p1[1]: sys.exit('error')
			N1 = self.pmap[e.p0]
			N2 = self.pmap[e.p1]
		else: #e.direct > 0
			if e.p1[1]>e.p0[1]: sys.exit('error')
			N1 = self.pmap[e.p1]
			N2 = self.pmap[e.p0]
		if e.base[0] < e.summit[0]: #x1<x2
			x1,y1 = e.base
			x2,y2 = e.summit
		else:
			x1,y1 = e.summit
			x2,y2 = e.base
		if x2==x1: print('Warning: Zero-Length Edges not eliminated clearly',e)
		a = 1.*(y2-y1)/(x2-x1)
		b = y1-a*x1
		x1 = round(x1,Voronoi.SLVround)
		if x1<0: x1 = 0
		elif x1==int(x1): x1 = int(x1)
		else: x1 = int(x1)+1
		x2=int(x2)+1
		if x2>self.ImgSizeX: x2 = self.ImgSizeX
		X = np.arange(x1,x2)
		Y = a*X+b
		for n in np.arange(X.size):
			Y[n] = round(Y[n],Voronoi.SLVround)
			if Y[n]<0: Y[n]=0
			nYn = int(Y[n])
			if Y[n] == nYn:
				self.pmap[X[n],Y[n]] = -1
				y1 = nYn-1
				y2 = nYn+1
			else:
				y1 = nYn
				y2 = nYn+1
			if y1>=0:
				N = self.pmap[X[n],y1]
				if N==0 or N>N1: self.pmap[X[n],y1] = N1
			if y2<self.ImgSizeY and self.pmap[X[n],y2]!=-1:
				N = self.pmap[X[n],y2]
				if N>=0 and N<N2: self.pmap[X[n],y2] = N2
		return

	def CalArea(self,**kwargs):
		self.ToCum2D=kwargs.pop('cum2d',self.ToCum2D)
		self.ToCalCtd = kwargs.get('calCentroid',self.ToCalCtd)
		self.ToCalTri=kwargs.pop('calTriangle',self.ToCalTri)
		self.ToRemoveEdgePoint = kwargs.pop('RemoveEdgePoint',self.ToRemoveEdgePoint)
		self.ToMarkEdgePoint = kwargs.pop('MarkEdgePoint',self.ToMarkEdgePoint)
		self.ToGetPoly = kwargs.pop('GetPoly',self.ToGetPoly)
		self.PolyFile = kwargs.pop('PolyFile',self.PolyFile)
		#for e in self.Edges.values():
		#	assert e.summit is not None
		print(color("Calculating Voronoi Cell Areas",32,0))
		if self.ToCalCtd: self.Wmap={(self.Px[n],self.Py[n]):np.zeros(2,dtype=np.float64) for n in np.arange(self.Px.size)}
		if type(self.pmap) is dict:
			self.Amap={(self.Px[n],self.Py[n]):np.float64(0) for n in np.arange(self.Px.size)}
			if self.ToGetPoly:
				self.Rmap={(self.Px[n],self.Py[n]):[] for n in np.arange(self.Px.size)} #Region map
				self.tang={(self.Px[n],self.Py[n]):[] for n in np.arange(self.Px.size)}
			P0 = np.array([e.p0 for e in self.Edges.values()])
			P1 = np.array([e.p1 for e in self.Edges.values()])
		else:
			self.Amap = np.zeros((self.ImgSizeX,self.ImgSizeY),dtype='float64')
			P0 = np.int32([e.p0 for e in self.Edges.values()])
			P1 = np.int32([e.p1 for e in self.Edges.values()])
		E0 = np.array([e.base   for e in self.Edges.values()]).round(Voronoi.SLVround)
		E1 = np.array([e.summit for e in self.Edges.values()]).round(Voronoi.SLVround)
		#KE = list(self.Edges.keys())
		#If keys, values and items views are iterated over with no intervening modifications to the dictionary, the order of items will directly correspond.
		for n in range(len(E0)):
			if (E0[n,0],E1[n,0]) in ((-0.5,self.RightLimit),(self.RightLimit,-0.5)) \
			or (E0[n,1],E1[n,1]) in ((-0.5,self.TopLimit),(self.TopLimit,-0.5)):
				print(P0[n],P1[n],E0[n],E1[n])

		DisPPList = lambda P0,P1: np.sqrt((P0[:,0]-P1[:,0])**2+(P0[:,1]-P1[:,1])**2)
		self.PPdis = DisPPList(P0,P1)

		#First calculate "Edge area list"
		#"Edge area": area of the triangle made of each edge and one of its according points
		length = DisPPList(E0,E1) 
		Earea = length * self.PPdis / 4.

		#Then calculate Parea
		#"Point area": area of each point (cell)
		for n in np.arange(P0.shape[0]):
			#self.Edges[KE[n]].length = length[n]
			#self.Edges[KE[n]].PPdist = self.PPdis[n]
			self.Amap[tuple(P0[n])] += Earea[n]
			self.Amap[tuple(P1[n])] += Earea[n]
			if self.ToGetPoly:
				at0=np.arctan2(E0[n,0]-P0[n,0],E0[n,1]-P0[n,1])
				at1=np.arctan2(E1[n,0]-P0[n,0],E1[n,1]-P0[n,1])
				atd=at0-at1
				if atd>np.pi: atd-=2*np.pi
				if atd<=-np.pi: atd+=2*np.pi
				if atd<=0:
					self.Rmap[tuple(P0[n])].append(E0[n])
					self.tang[tuple(P0[n])].append(at0)
				else:
					self.Rmap[tuple(P0[n])].append(E1[n])
					self.tang[tuple(P0[n])].append(at1)
				#if np.round(P0[n,0],3)==15.398 and np.round(P0[n,1],2)==608.84: print(P0[n],P1[n],E0[n],E1[n],at0,at1)
				at0=np.arctan2(E0[n,0]-P1[n,0],E0[n,1]-P1[n,1])
				at1=np.arctan2(E1[n,0]-P1[n,0],E1[n,1]-P1[n,1])
				atd=at0-at1
				if atd>np.pi: atd-=2*np.pi
				if atd<=-np.pi: atd+=2*np.pi
				if atd<=0:
					self.Rmap[tuple(P1[n])].append(E0[n])
					self.tang[tuple(P1[n])].append(at0)
				else:
					self.Rmap[tuple(P1[n])].append(E1[n])
					self.tang[tuple(P1[n])].append(at1)
				#if np.round(P1[n,0],3)==15.398 and np.round(P1[n,1],2)==608.84: print(P0[n],P1[n],E0[n],E1[n],at0,at1)
				'''
				at=np.arctan2(P1[n,0]-P0[n,0],P1[n,1]-P0[n,1])
				self.tang[tuple(P0[n])].append(at)
				if at>0: at-=np.pi
				else: at+=np.pi
				self.tang[tuple(P1[n])].append(at)
				if np.round(P0[n,0],3)==2876.136 and np.round(P0[n,1],2)==2762.47 or np.round(P1[n,0],3)==2876.136 and np.round(P1[n,1],2)==2762.47:
					print(P1[n],at)
				'''
			if self.ToCalCtd:
				self.Wmap[tuple(P0[n])] += Earea[n]*(E0[n]+E1[n])
				self.Wmap[tuple(P1[n])] += Earea[n]*(E0[n]+E1[n])

		EpL,EpR,EpB,EpT,VyL,VxB,VyR,VxT,XBeyond,XBelow,YBeyond,YBelow=self.FindBorderCells(P0,P1,E0,E1)
		if not self.ToRemoveEdgePoint and not self.ToMarkEdgePoint or self.ToGetPoly: 
			E0e,E1e,Pe,Ee=self.FillBorderCells(P0,P1,E0,E1,EpL,EpR,EpB,EpT,VyL,VxB,VyR,VxT,XBeyond,XBelow,YBeyond,YBelow)
			if self.ToGetPoly:
				for n in range(len(Pe)):
					at0=np.arctan2(E0e[n][0]-Pe[n][0],E0e[n][1]-Pe[n][1])
					at1=np.arctan2(E1e[n][0]-Pe[n][0],E1e[n][1]-Pe[n][1])
					atd=at0-at1
					if atd>np.pi: atd-=2*np.pi
					if atd<=-np.pi: atd+=2*np.pi
					if atd<=0:
						self.Rmap[tuple(Pe[n])].append(np.array(E0e[n]))
						self.tang[tuple(Pe[n])].append(at0)
					else:
						self.Rmap[tuple(Pe[n])].append(np.array(E1e[n]))
						self.tang[tuple(Pe[n])].append(at1)
					#if np.round(Pe[n][0],3)==15.398 and np.round(Pe[n][1],2)==608.84: print(Pe[n],at0,at1,atd)

		if self.ToGetPoly:
			for n in range(len(self.Px)):
				Px,Py=self.Px[n],self.Py[n]
				assert len(self.Rmap[(Px,Py)])==len(self.tang[(Px,Py)])
				self.Rmap[(Px,Py)] = np.array([self.Rmap[(Px,Py)][i] for i in np.argsort(self.tang[(Px,Py)])])-[self.OffSetX,self.OffSetY]
			if self.PolyFile is not None:
				with open(self.PolyFile,'w') as fout:
					print('image',file=fout)
					for n in range(len(self.Px)):
						print('polygon(',end='',file=fout)
						Px,Py=self.Px[n],self.Py[n]
						#if np.round(Px,3)==15.398 and np.round(Py,2)==608.84: print(self.Rmap[(Px,Py)])
						x,y=self.Rmap[(Px,Py)][0]
						print(f'{x:g},{y:g}',end='',file=fout)
						for x,y in self.Rmap[(Px,Py)][1:]: print(f',{x:g},{y:g}',end='',file=fout)
						print(f') # tag={{{len(self.Rmap[(Px,Py)])},{Px},{Py}}}',file=fout)

		if self.ToCum2D:
			#Calculate the 2D cummulative distribution
			#1. It's useless
			#2. Large cells at the image border introduce large error
			#3. Edge triangles not considered
			########################################
			#Divide the image into triangles (smallest pieces, number of triangles = number of edges)
			#Calculate the weight in each triangle
			Eweight0=np.zeros_like(Earea)
			Eweight1=np.zeros_like(Earea)
			for n in np.arange(P0.shape[0]):
				Eweight0[n] = Earea[n]/self.Amap[tuple(P0[n])]/len(self.Px)
				Eweight1[n] = Earea[n]/self.Amap[tuple(P1[n])]/len(self.Px)
			########################################
			fourcorner=np.array([[0,0],[0,self.ImgSizeY-1],[self.ImgSizeX-1,0],[self.ImgSizeX-1,self.ImgSizeY-1]])
			self.CumWmap={}
			for X,Y in np.vstack((E0,E1,fourcorner)): #all the cell vertexes and four image corners
				if (X,Y) in self.CumWmap.keys(): continue
				P0ok=(P0[:,0]<=X)&(P0[:,1]<=Y)
				P1ok=(P1[:,0]<=X)&(P1[:,1]<=Y)
				E0ok=(E0[:,0]<=X)&(E0[:,1]<=Y)
				E1ok=(E1[:,0]<=X)&(E1[:,1]<=Y)
				ok=np.vstack((P0ok,E0ok,E1ok,P1ok))
				#Sum all the triangles with at least two vertexes falling inside the integration box
				self.CumWmap[X,Y] = np.sum(Eweight0[np.sum(ok[:3],axis=0)>=2])+np.sum(Eweight1[np.sum(ok[1:],axis=0)>=2])
		################################################################################
		if self.ToCalTri:
			########################################
			#Divide the image into triangles (smallest pieces, number of triangles = number of edges)
			#Calculate the weight (number of points) in each triangle
			Eweight0=np.zeros_like(Earea)
			Eweight1=np.zeros_like(Earea)
			for n in np.arange(P0.shape[0]):
				Eweight0[n] = Earea[n]/self.Amap[tuple(P0[n])]
				Eweight1[n] = Earea[n]/self.Amap[tuple(P1[n])]
			with open(self.FileName+'_Triangles.dat','w') as fout:
				for n in np.arange(P0.shape[0]):
					print("%f %f %f %f %f %f %f %f" % (P0[n][0]-self.OffSetX,P0[n][1]-self.OffSetY,E0[n][0]-self.OffSetX,E0[n][1]-self.OffSetY,E1[n][0]-self.OffSetX,E1[n][1]-self.OffSetY,Earea[n],Eweight0[n]), file=fout)
					print("%f %f %f %f %f %f %f %f" % (P1[n][0]-self.OffSetX,P1[n][1]-self.OffSetY,E0[n][0]-self.OffSetX,E0[n][1]-self.OffSetY,E1[n][0]-self.OffSetX,E1[n][1]-self.OffSetY,Earea[n],Eweight1[n]), file=fout)
				if not self.ToRemoveEdgePoint:
					Eweighte=np.zeros_like(Ee)
					for n in np.arange(len(Pe)):
						Eweighte[n] = Ee[n]/self.Amap[tuple(Pe[n])]
					for n in np.arange(len(Pe)):
						print("%f %f %f %f %f %f %f %f" % (Pe[n][0]-self.OffSetX,Pe[n][1]-self.OffSetY,E0e[n][0]-self.OffSetX,E0e[n][1]-self.OffSetY,E1e[n][0]-self.OffSetX,E1e[n][1]-self.OffSetY,Ee[n],Eweighte[n]), file=fout)
				print('>> '+self.FileName+'_Triangles.dat')

		########################################################################
		if self.ToCalCtd:
			self.ctdvector=np.zeros((len(self.Px),2),dtype=np.float64)
			for n in range(len(self.Px)):
				Px,Py=self.Px[n],self.Py[n]
				self.Wmap[(Px,Py)] = (self.Wmap[(Px,Py)]/self.Amap[(Px,Py)]+[Px,Py])/3.
				self.ctdvector[n]= self.Wmap[(Px,Py)]-np.array([Px,Py])
			#print(list(self.Wmap.values())) print(list(self.Amap.keys()))
			#self.ctdvector=np.array(list(self.Wmap.values()))-np.array(list(self.Amap.keys()))
		########################################################################
		if self.ToMarkEdgePoint:
			for N in self.EdgePoint:
				self.Amap[self.Px[N-1],self.Py[N-1]] = -1
			print(color('EdgePoints marked as area=-1',34,1))
		if self.ToRemoveEdgePoint:
			for N in self.EdgePoint:
				self.Amap[self.Px[N-1],self.Py[N-1]] = 0
				self.Px[N-1]=-1
				self.Py[N-1]=-1
			self.Px=self.Px[self.Px>=0]
			self.Py=self.Py[self.Py>=0]
			print(color('EdgePoints removed',34,1))
	#END_OF_def CalArea(self):

	def FindBorderCells(self,P0,P1,E0,E1):
		################################################################################
		#Modify open cells at the image edge
		#Ignore the special cases of the cells at the four corners
		TriArea = lambda p1,p2,p3: abs(p2[0]*p3[1]+p3[0]*p1[1]+p1[0]*p2[1]-p1[0]*p3[1]-p2[0]*p1[1]-p3[0]*p2[1])/2.
		EE = np.vstack((E0,E1))

		########################################################################
		dd=0.5
		XBelow=-0.5-dd
		XBeyond=self.RightLimit+dd
		YBelow=-0.5-dd
		YBeyond=self.TopLimit+dd
		#L R B T [min,max]
		VyL = [YBeyond, YBelow]
		VxB = [XBeyond, XBelow]
		VyR = [YBeyond, YBelow]
		VxT = [XBeyond, XBelow]
		EpL = {}
		EpB = {}
		EpR = {}
		EpT = {}
		EdgePoint=[]
		################################################################################
		#Edge touching (crossing) two image edges (image corner)
		n = np.where((EE[:,0]==-0.5)&(EE[:,1]==-0.5))[0]	#bottom left
		if n.size>0:
			n = n[0]
			x1,y1 = P0[n-P0.shape[0]]
			x2,y2 = P1[n-P0.shape[0]]
			if x1 < x2:
				assert y1 > y2
				EpL[self.pmap[x1,y1]] = [-0.5]
				EpB[self.pmap[x2,y2]] = [-0.5]
			else:
				assert y1 < y2
				EpL[self.pmap[x2,y2]] = [-0.5]
				EpB[self.pmap[x1,y1]] = [-0.5]
			if Voronoi.debug: print(EE[n],'->',XBelow,YBelow)
			EE[n] = [XBelow,YBelow]
		n = np.where((EE[:,0]==-0.5)&(EE[:,1]==self.TopLimit))[0]	#top left
		if n.size>0:
			n = n[0]
			x1,y1 = P0[n-P0.shape[0]]
			x2,y2 = P1[n-P0.shape[0]]
			if x1 < x2:
				assert y1 < y2
				EpL[self.pmap[x1,y1]] = [self.TopLimit]
				EpT[self.pmap[x2,y2]] = [-0.5]
				#print self.pmap[x2,y2],EpT[self.pmap[x2,y2]]
			else:
				assert y1 > y2
				EpL[self.pmap[x2,y2]] = [self.TopLimit]
				EpT[self.pmap[x1,y1]] = [-0.5]
				#print self.pmap[x1,y1],EpT[self.pmap[x1,y1]]
			if Voronoi.debug: print(EE[n],'->',XBelow,YBeyond)
			EE[n] = [XBelow,YBeyond]
		n = np.where((EE[:,0]==self.RightLimit)&(EE[:,1]==-0.5))[0]	#bottom right
		if n.size>0:
			n = n[0]
			x1,y1 = P0[n-P0.shape[0]]
			x2,y2 = P1[n-P0.shape[0]]
			if x1 < x2:
				assert y1 < y2
				EpB[self.pmap[x1,y1]] = [self.RightLimit]
				EpR[self.pmap[x2,y2]] = [-0.5]
			else:
				assert y1 > y2
				EpB[self.pmap[x2,y2]] = [self.RightLimit]
				EpR[self.pmap[x1,y1]] = [-0.5]
			if Voronoi.debug: print(EE[n],'->',XBeyond,YBelow)
			EE[n] = [XBeyond,YBelow]
		n = np.where((EE[:,0]==self.RightLimit)&(EE[:,1]==self.TopLimit))[0]	#top right
		if n.size>0:
			n = n[0]
			x1,y1 = P0[n-P0.shape[0]]
			x2,y2 = P1[n-P0.shape[0]]
			if x1 < x2:
				assert y1 > y2
				EpT[self.pmap[x1,y1]] = [self.RightLimit]
				EpR[self.pmap[x2,y2]] = [self.TopLimit]
				#print self.pmap[x1,y1],EpT[self.pmap[x1,y1]]
			else:
				assert y1 < y2
				EpT[self.pmap[x2,y2]] = [self.RightLimit]
				EpR[self.pmap[x1,y1]] = [self.TopLimit]
				#print self.pmap[x2,y2],EpT[self.pmap[x2,y2]]
			if Voronoi.debug: print(EE[n],'->',XBeyond,YBeyond)
			EE[n] = [XBeyond,YBeyond]
		################################################################################
		#Edges touching (crossing) one image edge
		if Voronoi.debug: EPfile = open('EdgePoint.reg','w')
		for n in np.where(EE[:,0]==-0.5)[0]:	#left EpL
			#print('EpL',EE[n])
			N = self.pmap[tuple(P0[n-P0.shape[0]])]
			EdgePoint.append(N)
			EpL[N] = EpL.get(N,[])+[EE[n,1]]
			N = self.pmap[tuple(P1[n-P0.shape[0]])]
			EdgePoint.append(N)
			EpL[N] = EpL.get(N,[])+[EE[n,1]]
			VyL[0] = min(VyL[0],EE[n,1])
			VyL[1] = max(VyL[1],EE[n,1])
		for n in np.where(EE[:,1]==-0.5)[0]:	#bottom EpB
			#print(n,EE[n], EE[n-P0.shape[0]], E0[n],E1[n], P0[n-P0.shape[0]], P1[n-P0.shape[0]])
			#print(n,EE[n], EE[n-P0.shape[0]], E0[n-P0.shape[0]],E1[n-P0.shape[0]], P0[n-P0.shape[0]], P1[n-P0.shape[0]])
			N = self.pmap[tuple(P0[n-P0.shape[0]])]
			EdgePoint.append(N)
			#print(N,EpB.get(N,[]),end='->')
			EpB[N] = EpB.get(N,[])+[EE[n,0]]
			#print(EpB[N])
			N = self.pmap[tuple(P1[n-P0.shape[0]])]
			EdgePoint.append(N)
			#print(N,EpB.get(N,[]),end='->')
			EpB[N] = EpB.get(N,[])+[EE[n,0]]
			#print(EpB[N])
			VxB[0] = min(VxB[0],EE[n,0])
			VxB[1] = max(VxB[1],EE[n,0])
		for n in np.where(EE[:,0]==self.RightLimit)[0]:	#right EpR
			N = self.pmap[tuple(P0[n-P0.shape[0]])]
			EdgePoint.append(N)
			EpR[N] = EpR.get(N,[])+[EE[n,1]]
			#print('a',N,EE[n],self.Px[N-1],self.Py[N-1])
			if Voronoi.debug: print('line(%.2f,%.2f,%.2f,%.2f)' %(self.Py[N-1]+1,self.Px[N-1]+1,EE[n,1]+1,EE[n,0]+1), file=EPfile)
			N = self.pmap[tuple(P1[n-P0.shape[0]])]
			EdgePoint.append(N)
			EpR[N] = EpR.get(N,[])+[EE[n,1]]
			#print('bb',N,EE[n],self.Px[N-1],self.Py[N-1])
			#if N==3843: print( P0[n-len(P0)],P1[n-len(P0)],list(self.Edges.values())[n-len(P0)]) #!!!
			if Voronoi.debug: print('line(%.2f,%.2f,%.2f,%.2f)' %(self.Py[N-1]+1,self.Px[N-1]+1,EE[n,1]+1,EE[n,0]+1), file=EPfile)
			VyR[0] = min(VyR[0],EE[n,1])
			VyR[1] = max(VyR[1],EE[n,1])
		for n in np.where(EE[:,1]==self.TopLimit)[0]:	#top EpT
			N = self.pmap[tuple(P0[n-P0.shape[0]])]
			EdgePoint.append(N)
			EpT[N] = EpT.get(N,[])+[EE[n,0]]
			N = self.pmap[tuple(P1[n-P0.shape[0]])]
			EdgePoint.append(N)
			EpT[N] = EpT.get(N,[])+[EE[n,0]]
			VxT[0] = min(VxT[0],EE[n,0])
			VxT[1] = max(VxT[1],EE[n,0])
		self.EdgePoint = np.unique(EdgePoint)
		if Voronoi.debug:
			for N in self.EdgePoint:
				print("circle(%d,%d,1) # color=red" % (self.Py[N-1]+1-self.OffSetY,self.Px[N-1]+1-self.OffSetX), file=EPfile)
			EPfile.close()
			print(">> "+EPfile.name)
		return EpL,EpR,EpB,EpT,VyL,VxB,VyR,VxT,XBeyond,XBelow,YBeyond,YBelow
	#END_OF_def FindBorderCells(self,P0,P1,E0,E1):

	def FillBorderCells(self,P0,P1,E0,E1,EpL,EpR,EpB,EpT,VyL,VxB,VyR,VxT,XBeyond,XBelow,YBeyond,YBelow):
		#print('EpL',EpL)
		#print('EpR',EpR)
		#print('EpB',EpB)
		#print('EpT',EpT)
		#print('VyL',VyL)
		#print('VyR',VyR)
		#print('VxT',VxT)
		#print('VxB',VxB)
		if len(EpR)==len(EpT)==0:
			assert VyR == [YBeyond,YBelow] and VxT == [XBeyond,XBelow]
			n = np.argmax(self.Px+self.Py)
			N = self.pmap[self.Px[n],self.Py[n]]
			EpR[N] = [-0.5,self.TopLimit] #lower point -0.5 to be filled below with EpB. NO. don't
			EpT[N] = [-0.5,self.RightLimit] #left point -0.5 to be filled below with EpL. NO. don't
			assert len(EpL[N])==1
			EpL[N].append(self.TopLimit)
			assert len(EpB[N])==1
			EpB[N].append(self.RightLimit)
		elif len(EpL)==len(EpT)==0:
			assert VyL == [YBeyond,YBelow] and VxT == [XBeyond,XBelow]
			n = np.argmax(-self.Px+self.Py)
			N = self.pmap[self.Px[n],self.Py[n]]
			EpL[N] = [-0.5,self.TopLimit]
			EpT[N] = [-0.5,self.RightLimit]
			assert len(EpR[N])==1
			EpR[N].append(self.TopLimit)
			assert len(EpB[N])==1
			EpB[N].append(-0.5)
		elif len(EpR)==len(EpB)==0:
			assert VyR == [YBeyond,YBelow] and VxB == [XBeyond,XBelow]
			n = np.argmax(self.Px-self.Py)
			N = self.pmap[self.Px[n],self.Py[n]]
			EpR[N] = [-0.5,self.TopLimit]
			EpB[N] = [-0.5,self.RightLimit]
			assert len(EpL[N])==1
			EpL[N].append(-0.5)
			assert len(EpT[N])==1
			EpT[N].append(self.RightLimit)
		elif len(EpL)==len(EpB)==0:
			assert VyL == [YBeyond,YBelow] and VxB == [XBeyond,XBelow]
			n = np.argmax(-self.Px-self.Py)
			N = self.pmap[self.Px[n],self.Py[n]]
			EpL[N] = [-0.5,self.TopLimit]
			EpB[N] = [-0.5,self.RightLimit]
			assert len(EpR[N])==1
			EpR[N].append(-0.5)
			assert len(EpT[N])==1
			EpT[N].append(-0.5)
		else:
			# For such case where there is no crossing point in one edge of the image (corner quadrangle cell)
			if len(EpL)==0:
				assert VyL == [YBeyond,YBelow]
				n = np.argmin(self.Px)
				N = self.pmap[self.Px[n],self.Py[n]]
				assert len(EpB[N])==1 and len(EpT[N])==1 and EpB[N][0]==VxB[0] and EpT[N][0]==VxT[0]
				EpL[N] = [-0.5, self.TopLimit]
				EpB[N].append(-0.5)
				EpT[N].append(-0.5)
			if len(EpB)==0:
				assert VxB == [XBeyond,XBelow]
				n = np.argmin(self.Py)
				N = self.pmap[self.Px[n],self.Py[n]]
				#print(EpL.get(N,None),EpR.get(N,None),EpL.get(N,None),'0=',VyL,EpR.get(N,None),'0=',VyR)
				assert len(EpL[N])==1 and len(EpR[N])==1 and EpL[N][0]==VyL[0] and EpR[N][0]==VyR[0]
				EpB[N] = [-0.5, self.RightLimit]
				EpL[N].append(-0.5)
				EpR[N].append(-0.5)
			if len(EpR)==0:
				assert VyR == [YBeyond,YBelow]
				n = np.argmax(self.Px)
				N = self.pmap[self.Px[n],self.Py[n]]
				assert len(EpB[N])==1 and len(EpT[N])==1 and EpB[N][0]==VxB[1] and EpT[N][0]==VxT[1]
				EpR[N] = [-0.5, self.TopLimit]
				EpB[N].append(self.RightLimit)
				EpT[N].append(self.RightLimit)
			if len(EpT)==0:
				assert VxT == [XBeyond,XBelow]
				n = np.argmax(self.Py)
				N = self.pmap[self.Px[n],self.Py[n]]
				assert len(EpL[N])==1 and len(EpR[N])==1 and EpL[N][0]==VyL[1] and EpR[N][0]==VyR[1]
				EpT[N] = [-0.5, self.RightLimit]
				EpL[N].append(self.TopLimit)
				EpR[N].append(self.TopLimit)
		########################################################################
		#corner triangle cell and quadrangle cell
		#make sure each Point has two Edge points
		#print('EpL',EpL)
		#print('EpR',EpR)
		#print('EpB',EpB)
		#print('EpT',EpT)
		#print('VyL',VyL)
		#print('VyR',VyR)
		#print('VxT',VxT)
		#print('VxB',VxB)
		for N in EpL:
			x,y = self.Px[N-1],self.Py[N-1]
			if len(EpL[N])==1:
				if EpL[N][0]==-0.5:
					assert [-0.5] in EpB.values()
					if VyL[0]<=self.TopLimit:
						EpL[N].append(VyL[0])
					else:
						EpL[N].append(self.TopLimit)
						assert len(EpT[N])==1
						EpT[N].append(-0.5)
				elif EpL[N][0]==self.TopLimit:
					assert [-0.5] in EpT.values()
					if VyL[1]>=-0.5:
						EpL[N].append(VyL[1])
					else:
						EpL[N].append(-0.5)
						assert len(EpB[N])==1
						EpB[N].append(-0.5)
				#elif EpL[N][0]==VyL[0] and N in EpB.keys(): #divide into 3 cases below
				#	#N refers to the left bottom pixel
				#	#print(x,y,N,EpL[N],EpB[N],VxB)
				#	assert len(EpB[N])==1 and EpB[N][0]==min(VxB[0],self.RightLimit)
				#	EpL[N].append(-0.5)
				#	EpB[N].append(-0.5)
				elif EpL[N][0]==VyL[0] and N in EpB.keys() and N in EpT.keys() and N in EpL.keys():
					#Voronoi_second_corner_pixel.png, touching 3 edges but not the very-corner pixel
					assert len(EpB[N])==2
					assert len(EpT[N])==1 and EpT[N][0]==min(VxT[0],self.RightLimit)
					assert len(EpL[N])==1 and EpL[N][0]==max(VyL[-1],-0.5)
					EpL[N].append(self.TopLimit)
					EpT[N].append(-0.5)
				elif EpL[N][0]==VyL[0] and N in EpB.keys() and N in EpL.keys() and len(EpB[N])==len(EpL[N])==2:
					#like Voronoi_second_corner_pixel.png, but touching only B and L and not T
					pass
				elif EpL[N][0]==VyL[0] and N in EpB.keys():
					assert len(EpB[N])==1 and EpB[N][0]==min(VxB[0],self.RightLimit)
					EpL[N].append(-0.5)
					EpB[N].append(-0.5)
				elif EpL[N][0]==VyL[1] and N in EpT.keys():
					#Only one cross point in this edge. N refers to the left top pixel
					#print(N,EpT[N],EpL[N],EpT[N])
					#print(VxT[0],self.RightLimit)
					assert len(EpT[N])==1 and EpT[N][0]==min(VxT[0],self.RightLimit)
					EpL[N].append(self.TopLimit)
					EpT[N].append(-0.5)
				else:
					print('EpL',x,y,EpL[N])
					print('ERROR--------------------------')
					raise RuntimeError('ERROR: Edge point Left')
		for N in EpR:
			x,y = self.Px[N-1],self.Py[N-1]
			if len(EpR[N])==1:
				if EpR[N][0]==-0.5:
					assert [self.RightLimit] in EpB.values()
					if VyR[0]<=self.TopLimit:
						EpR[N].append(VyR[0])
					else:
						EpR[N].append(self.TopLimit)
						assert len(EpT[N])==1
						EpT[N].append(self.RightLimit)
				elif EpR[N][0]==self.TopLimit:
					assert [self.RightLimit] in EpT.values()
					if VyR[1]>=-0.5:
						EpR[N].append(VyR[1])
					else:
						EpR[N].append(-0.5)
						assert len(EpB[N])==1
						EpT[N].append(self.RightLimit)
				#elif EpR[N][0]==VyR[0] and N in EpB.keys():
				#	assert len(EpB[N])==1 and EpB[N][0]==max(VxB[1],-0.5)
				#	EpR[N].append(-0.5)
				#	EpB[N].append(self.RightLimit)
				elif EpR[N][0]==VyR[0] and N in EpB.keys() and N in EpR.keys() and N in EpT.keys():
					#Voronoi_second_corner_pixel.png, touching 3 edges but not the very-corner pixel
					assert len(EpB[N])==2
					assert len(EpT[N])==1 and EpT[N][0]==max(VxT[-1],-0.5)
					assert len(EpR[N])==1 and EpR[N][0]==max(VyR[-1],-0.5)
					EpR[N].append(self.TopLimit)
					EpT[N].append(self.RightLimit)
				elif EpR[N][0]==VyR[0] and N in EpB.keys() and N in EpR.keys() and len(EpB[N])==len(EpR[N])==2:
					#like Voronoi_second_corner_pixel.png, but touching only B and L and not T
					pass
				elif EpR[N][0]==VyR[0] and N in EpB.keys():
					assert len(EpB[N])==1 and EpB[N][0]==max(VxB[1],-0.5)
					EpR[N].append(-0.5)
					EpB[N].append(self.RightLimit)
				elif EpR[N][0]==VyR[1] and N in EpT.keys():
					assert len(EpT[N])==1 and EpT[N][0]==max(VxT[1],-0.5)
					EpR[N].append(self.TopLimit)
					EpT[N].append(self.RightLimit)
				else:
					print('EpR',x,y,N,EpR[N])
					#'''
					n= np.where((P0[:,0]==x)&(P0[:,1]==y))
					if len(n)>0:
						n=n[0][0]
						print('P0',n,P0[n]) #,EE[n],EE[n-P0.shape[0]])
					n= np.where((P1[:,0]==x)&(P1[:,1]==y))
					if len(n)>0:
						n=n[0][0]
						print('P1',n,P1[n]) #,EE[n],EE[n-P1.shape[0]])
					#'''
					raise RuntimeError('ERROR: Edge point Right')
		########################################################################
		#Save edge triangles in E0e,E1e,Pe,Ee. Add their areas to self.Amap
		E0e=[] #vertex at the border
		E1e=[] #vertex at the border
		Pe=[] #The EdgePoint of each EdgeTriangle.
		#The number can be larger than that of self.EdgePoint, because one EdgePoint can be shared by more than one triangles.
		Ee=[] #area of each edge triangle
		for N in EpL:
			x,y = self.Px[N-1],self.Py[N-1]
			#if len(EpL[N])!=2: print y+1,x+1,EpL[N]
			assert len(EpL[N])==2
			t_area = (x+0.5)*abs(EpL[N][0]-EpL[N][1])/2.
			#print('Left',x,y,EpL[N])
			#if N in EpB: print('	bottom',EpB[N])
			#if N in EpT: print('	top',EpT[N])
			self.Amap[x,y] += t_area
			if self.ToCalCtd:
				self.Wmap[x,y] += t_area*np.array([-0.5*2,EpL[N][0]+EpL[N][1]]) #edge_left: x=minimum, two y values
			E0e.append([-0.5,EpL[N][0]])
			E1e.append([-0.5,EpL[N][1]])
			Pe.append([x,y])
			Ee.append((x+0.5)*abs(EpL[N][0]-EpL[N][1])/2.)
			#print "circle(%.1f,%.1f,1.3) # text={%.3f}" %(y+1,x+1,(x+0.5)*abs(EpL[N][0]-EpL[N][1])/2.)
		for N in EpB:
			x,y = self.Px[N-1],self.Py[N-1]
			#if len(EpB[N])!=2: print y+1,x+1,EpB[N]
			assert len(EpB[N])==2
			t_area = (y+0.5)*abs(EpB[N][0]-EpB[N][1])/2.
			#print('Bottom',x,y,EpB[N])
			#if N in EpL: print('	left',EpL[N])
			#if N in EpR: print('	right',EpR[N])
			self.Amap[x,y] += t_area
			if self.ToCalCtd:
				self.Wmap[x,y] += t_area*np.array([EpB[N][0]+EpB[N][1],-0.5*2]) #edge_bottom: two x values, y=minimum
			#print(x,y,t_area)
			E0e.append([EpB[N][0],-0.5])
			E1e.append([EpB[N][1],-0.5])
			Pe.append([x,y])
			Ee.append((y+0.5)*abs(EpB[N][0]-EpB[N][1])/2.)
			#print "circle(%.1f,%.1f,1.3) # text={%.3f}" %(y+1,x+1,(y+0.5)*abs(EpB[N][0]-EpB[N][1])/2.)
		for N in EpR:
			x,y = self.Px[N-1],self.Py[N-1]
			#if len(EpR[N])!=2: print y+1,x+1,EpR[N]
			assert len(EpR[N])==2
			t_area = (self.RightLimit-x)*abs(EpR[N][0]-EpR[N][1])/2.
			self.Amap[x,y] += t_area
			#print('Right',x,y,EpR[N])
			#if N in EpB: print('	bottom',EpB[N])
			#if N in EpT: print('	top',EpT[N])
			if self.ToCalCtd:
				self.Wmap[x,y] += t_area*np.array([(self.RightLimit)*2,EpR[N][0]+EpR[N][1]]) #edge_right: x=maximum, two y values
			E0e.append([self.RightLimit,EpR[N][0]])
			E1e.append([self.RightLimit,EpR[N][1]])
			Pe.append([x,y])
			Ee.append((self.RightLimit-x)*abs(EpR[N][0]-EpR[N][1])/2.)
			#print "circle(%.1f,%.1f,1.3) # text={%.3f}" %(y+1,x+1,(self.RightLimit-x)*abs(EpR[N][0]-EpR[N][1])/2.)
		for N in EpT:
			x,y = self.Px[N-1],self.Py[N-1]
			#if len(EpT[N])!=2: print y+1,x+1,EpT[N],N
			assert len(EpT[N])==2
			t_area = (self.TopLimit-y)*abs(EpT[N][0]-EpT[N][1])/2.
			#print('Top',x,y,EpT[N])
			#if N in EpL: print('	left',EpL[N])
			#if N in EpR: print('	right',EpR[N])
			self.Amap[x,y] += t_area
			if self.ToCalCtd:
				self.Wmap[x,y] += t_area*np.array([EpT[N][0]+EpT[N][1],(self.TopLimit)*2]) #edge_top, two x values, y=maximum
			E0e.append([EpT[N][0],self.TopLimit])
			E1e.append([EpT[N][1],self.TopLimit])
			Pe.append([x,y])
			Ee.append((self.TopLimit-y)*abs(EpT[N][0]-EpT[N][1])/2.)
			#print "circle(%.1f,%.1f,1.3) # text={%.3f}" %(y+1,x+1,(self.TopLimit-y)*abs(EpT[N][0]-EpT[N][1])/2.)
		if Voronoi.debug:
			with open('EdgeTriangle.reg','w') as freg:
				for n in range(len(Pe)):
					print("polygon(%f,%f,%f,%f,%f,%f) # tag={%f}" % (E0e[n][1]+1-self.OffSetY,E0e[n][0]+1-self.OffSetX,E1e[n][1]+1-self.OffSetY,E1e[n][0]+1-self.OffSetX,Pe[n][1]+1-self.OffSetY,Pe[n][0]+1-self.OffSetX,Ee[n]), file=freg)
				print(">> "+freg.name)
		return E0e,E1e,Pe,Ee
	#END_OF_def FillBorderCells(self,P0,P1,E0,E1,EpL,EpR,EpB,EpT,VyL,VxB,VyR,VxT):
	@staticmethod
	def anscombe(x):
		assert np.all(x>=0)
		y = 2.*np.sqrt(x+3./8.)
		return y
	@staticmethod
	def inverse_anscombe(z):
		'''
		Reference: Makitalo, M., & Foi, A. (2011). A closed-form approximation of the exact unbiased inverse of the Anscombe variance-stabilizing transformation. Image Processing.
		'''
		return 0.25*np.power(z,2) + 0.25*np.sqrt(1.5)*np.power(z,-1) - 11./8.*np.power(z,-2) + 5./8.*np.sqrt(1.5)*np.power(z,-3) -0.125

	@staticmethod
	def PixelOutflow(f1map,Nmap,P0,P1,SmoothFactor):
		if f1map is dict:
			Smap={k:SmoothFactor*f1map[k] for k in f1map.keys()}
			f2map={k:(1-SmoothFactor)*f1map[k]/Nmap[k] for k in f1map.keys()}
		else:
			Smap=SmoothFactor*f1map
			f2map=(1-SmoothFactor)*f1map
			haveP=Nmap>0
			f2map[haveP]/=Nmap[haveP]
		for n in np.arange(P0.shape[0]):
			Smap[tuple(P0[n])] += f2map[tuple(P1[n])]
			Smap[tuple(P1[n])] += f2map[tuple(P0[n])]
		return Smap
	@staticmethod
	def PixelInflow(dmap,Nmap,P0,P1,SmoothFactor):
		if dmap is dict:
			Smap={k:np.float64(0) for k in dmap.keys()}
		else:
			Smap=np.zeros_like(dmap)
			Vmap=np.zeros_like(dmap)
		for n in np.arange(P0.shape[0]):
			Smap[tuple(P0[n])] += dmap[tuple(P1[n])]
			Smap[tuple(P1[n])] += dmap[tuple(P0[n])]
		if dmap is dict:
			for k in Smap: Smap[k]=dmap[ok]*SmoothFactor+Smap[k]*(1-SmoothFactor)/Nmap[k]
		else:
			haveP=Nmap>0
			Smap[haveP]/=Nmap[haveP]
			if 1:
				for n in np.arange(P0.shape[0]):
					Vmap[tuple(P0[n])] += (dmap[tuple(P1[n])]-Smap[tuple(P1[n])])**2
					Vmap[tuple(P1[n])] += (dmap[tuple(P0[n])]-Smap[tuple(P0[n])])**2
				Vmap[haveP]=np.sqrt(Vmap[haveP]/(Nmap[haveP]-1))/Smap[haveP]
			Smap*=(1-SmoothFactor)
			Smap+=dmap*SmoothFactor
		if 1: return Smap,Vmap
		else: return Smap
	def SplitImage(self,img,MaxBkgDensity,MaxBkgCountPP):
		img1=img.copy()
		maximg=np.int32(np.ceil(self.Amap*MaxBkgDensity))
		high=img>maximg
		img1[high]=maximg[high]
		img1[img1>MaxBkgCountPP]=MaxBkgCountPP
		assert np.all(np.argwhere(img>0)==np.argwhere(img1>0))
		assert np.all(img1>=0)
		return img1,img-img1
	def SmoothBackground(self,Cmap,**kwargs):
		MaxBkgDensity=kwargs.get('MaxBkgDensity',1.) #MaxBkgDensity=1 is bad for deep fields like SEP
		MaxBkgCountPP=kwargs.get('MaxBkgCountPP',3)
		assert np.all(self.Amap[self.Px,self.Py]>0) #self.CalArea done
		assert np.all(Cmap[self.Px,self.Py]>0)
		assert type(self.pmap) is not dict
		img1,img2=self.SplitImage(Cmap,MaxBkgDensity,MaxBkgCountPP)
		self.SmoothDensity(Cmap=img1,**kwargs)
		dmap=self.Smap.copy() #full densiy map
		haveP=img2>0
		dmap[haveP]+=(img2[haveP]/self.Amap[haveP])
		return dmap,self.Smap
	def SmoothDensity(self,Cmap=None,**kwargs):
		'''
		1. Smooth density map, not count map. Although the latter keeps the total counts, it is meaningless in the sense of Voronoi smoothing.
		2. Outflowing from each pixel results in a non-uniform map proportional to Nmap. Both Outflow and Inflow lead to larger total counts.
		If Cmap=None, 1 is adopted, because the Voronoi class keeps only point positions but not point values.
		'''
		print(color("Voronoi Smooothing",32,0))
		assert np.all(self.Amap[self.Px,self.Py]>0) #self.CalArea done
		assert np.all(Cmap[self.Px,self.Py]>0)
		self.SmoothNumber = kwargs.pop('SmoothNumber',self.SmoothNumber)
		self.SmoothFactor = kwargs.get('SmoothFactor',self.SmoothFactor)
		self.ToRemoveEdgePoint = kwargs.pop('RemoveEdgePoint',self.ToRemoveEdgePoint)
		self.ToMarkEdgePoint = kwargs.pop('MarkEdgePoint',self.ToMarkEdgePoint)
		Mask=kwargs.pop('Mask',None)
		CutMask=kwargs.pop('CutMask',False)
		if type(self.pmap) is dict:
			P0 = np.array([e.p0 for e in self.Edges.values()])
			P1 = np.array([e.p1 for e in self.Edges.values()])
			if Cmap is None: Cmap={k:1. for k in self.Amap.keys()}
			assert all(Cmap[k]>0 for k in self.Amap.values())
			Nmap={k:0 for k in self.Amap.keys()}
			dmap={k:Cmap[k]/self.Amap[k] for k in self.Amap.keys()}
		else:
			P0 = np.int32([e.p0 for e in self.Edges.values()])
			P1 = np.int32([e.p1 for e in self.Edges.values()])
			if Cmap is None:
				Cmap = np.zeros_like(self.Amap)
				Cmap[self.Px,self.Py] = 1
			Nmap=np.zeros(Cmap.shape,dtype=np.int16)
			dmap=np.divide(Cmap,self.Amap,out=np.zeros_like(self.Amap),where=self.Amap>0)
		if Mask is not None:
			assert type(Mask)==type(self.Amap)
			good = (Mask[P0[:,0],P0[:,1]]>0)&(Mask[P1[:,0],P1[:,1]]>0)
			bad  = (Mask[P0[:,0],P0[:,1]]<=0)&(Mask[P1[:,0],P1[:,1]]<=0)
			P0m=P0[bad]
			P1m=P1[bad]
			if CutMask:
				P0 = P0[good]
				P1 = P1[good]
		for n in np.arange(P0.shape[0]):
			Nmap[tuple(P0[n])]+=1
			Nmap[tuple(P1[n])]+=1
		ANS=True #It only slightly make the result flatter (~1%)
		if ANS: dmap=self.anscombe(dmap)
		#fits.writeto('dmap.fits',dmap)
		if 1:
			#self.Smap,self.Vmap=self.PixelInflow(dmap,Nmap,P0,P1,self.SmoothFactor)
			self.Smap=dmap
			for i in range(self.SmoothNumber):
				self.Smap,self.Vmap=self.PixelInflow(self.Smap,Nmap,P0,P1,self.SmoothFactor)
				#fits.writeto(f'Smap{i}.fits',griddata(np.vstack((self.Px,self.Py)).T,self.Smap[self.Px,self.Py],tuple(np.indices((self.ImgSizeX,self.ImgSizeY))),method='linear',fill_value=np.median(self.Smap[self.Px,self.Py])))
				#fits.writeto(f'Vmap{i}.fits',griddata(np.vstack((self.Px,self.Py)).T,self.Vmap[self.Px,self.Py],tuple(np.indices((self.ImgSizeX,self.ImgSizeY))),method='linear',fill_value=np.median(self.Vmap[self.Px,self.Py])))

				#print(i,np.sum(Cmap),np.sum(self.Smap*self.Amap)) #the total counts becomes larger and larger
		else:
			self.Smap=self.PixelInflow(dmap,Nmap,P0,P1,self.SmoothFactor)
			for i in range(1,self.SmoothNumber):
				self.Smap=self.PixelInflow(self.Smap,Nmap,P0,P1,self.SmoothFactor)
		if ANS: self.Smap=self.inverse_anscombe(self.Smap)
		#fits.writeto('Smap3.fits',self.Smap)
		scmap=self.Smap*self.Amap
		sc=np.log10(scmap[scmap>0])
		for i in range(3):
			upper3sigma=np.std(sc)*3+np.median(sc)
			sc=sc[sc<upper3sigma]
		if Mask is not None: reg=(Cmap>0)&(np.log10(scmap)<upper3sigma)&(Mask>0)
		else: reg=(Cmap>0)&(np.log10(scmap)<upper3sigma)
		f=np.sum(Cmap[reg])/np.sum(self.Smap[reg]*self.Amap[reg])
		self.Smap*=f
		self.Smap[Cmap<=0]=0
		#fits.writeto('Smap4.fits',self.Smap)
		if Mask is not None:
			self.Smap[P0m[0],P0m[1]]=dmap[P0m[0],P0m[1]]
			self.Smap[P1m[0],P1m[1]]=dmap[P1m[0],P1m[1]]
		if self.ToRemoveEdgePoint:
			for N in self.EdgePoint:
				self.Amap[self.Px[N-1],self.Py[N-1]] = 0
				self.Smap[self.Px[N-1],self.Py[N-1]] = 0
				self.Px[N-1]=-1
				self.Py[N-1]=-1
			self.Px=self.Px[self.Px>=0]
			self.Py=self.Py[self.Py>=0]
			print(color('EdgePoints removed',34,1))

	def CalDVD(self,**kwargs):
		#Still has a problem as shown in problem_DVD_20230502.png !!!!
		#pixel at the vertex is set by later (larger number) edges although the earlier one is correct
		ImgSizeX=self.ImgSizeX
		ImgSizeY=self.ImgSizeY
		#print(ImgSizeX,ImgSizeY)
		NP=2 #2 pixels padding
		fig = pl.Figure(figsize=((ImgSizeX+NP*2)/100, (ImgSizeY+NP*2)/100), dpi=100)
		canvas = FigureCanvas(fig)
		ax = fig.add_axes([0,0,1,1])
		for n in tqdm(np.arange(0,len(self.Px))):
			polygon=np.array(self.Rmap[(self.Px[n],self.Py[n])])+NP+[self.OffSetX,self.OffSetY]
			ax.fill(polygon[:,0],polygon[:,1],color=INT2RGB(n+1),edgecolor='None',aa=False)
			#it fills each polygon correctly, but some pixels whose center is just outside the border are included
			#so the cells filled later will be larger and the ones done earlier will be smaller
		ax.set_xlim(-0.5,ImgSizeX+-0.5+NP*2)
		ax.set_ylim(-0.5,ImgSizeY+-0.5+NP*2)
		canvas.draw()
		assert canvas.get_width_height()==(fig.figure.get_window_extent().width,fig.figure.get_window_extent().height)
		d=np.int32(np.frombuffer(canvas.tostring_rgb(), dtype='uint8').reshape(ImgSizeY+NP*2,ImgSizeX+NP*2,3))
		i=imgRGB2INT(d)
		self.filledmap=i[ImgSizeY+NP-1:NP-1:-1,NP:ImgSizeX+NP]
		#fits.writeto('tmp1.fits',self.filledmap,overwrite=True)

		ax.clear()
		#st=time.time()
		for n in tqdm(self.Edges.keys()):
			e=self.Edges[n] #positions in Python format
			ax.plot((e.base[0]+NP,e.summit[0]+NP),(e.base[1]+NP,e.summit[1]+NP),color=INT2RGB(n),lw=0.1,aa=False)
		#It's funny that the loop above is faster than the no-loop method below.
		#I guess the loop in the plot() function is inefficient
		#lineset=[((self.Edges[n].base[0]+NP,self.Edges[n].summit[0]+NP),(self.Edges[n].base[1]+NP,self.Edges[n].summit[1]+NP),INT2HEX(n)) for n in self.Edges.keys()]
		#ax.plot(*tqdm(itertools.chain.from_iterable(lineset)),lw=0.1,aa=False)
		#print(time.time()-st)
		ax.set_xlim(-0.5,ImgSizeX+-0.5+NP*2)
		ax.set_ylim(-0.5,ImgSizeY+-0.5+NP*2)
		canvas.draw()
		assert canvas.get_width_height()==(fig.figure.get_window_extent().width,fig.figure.get_window_extent().height)
		d=np.int32(np.frombuffer(canvas.tostring_rgb(), dtype='uint8').reshape(ImgSizeY+NP*2,ImgSizeX+NP*2,3))
		i=imgRGB2INT(d)
		lmap=i[ImgSizeY+NP-1:NP-1:-1,NP:ImgSizeX+NP]
		#fits.writeto('tmp2.fits',lmap,overwrite=True)
	
		if 0: #line drawn later will overwrite the correct point. try to identify such pixels in order to make correction.
			#But this plotting does not yet work. It's unclear which pixel is plotted
			ax.clear()
			pset=[([e.base[0]+NP],[e.base[1]+NP],'b,',[e.summit[0]+NP],[e.summit[1]+NP],'r,') for e in self.Edges.values()] #seems "," pixel marker often offset
			#pset=[([e.base[0]+NP],[e.base[1]+NP],'bs',[e.summit[0]+NP],[e.summit[1]+NP],'rs') for e in self.Edges.values()] #seems "," pixel marker often offset
			ax.plot(*tqdm(itertools.chain.from_iterable(pset)),aa=False,ms=0.1)
			ax.set_xlim(-0.5,ImgSizeX+-0.5+NP*2)
			ax.set_ylim(-0.5,ImgSizeY+-0.5+NP*2)
			canvas.draw()
			assert canvas.get_width_height()==(fig.figure.get_window_extent().width,fig.figure.get_window_extent().height)
			d=np.int32(np.frombuffer(canvas.tostring_rgb(), dtype='uint8').reshape(ImgSizeY+NP*2,ImgSizeX+NP*2,3))
			i=imgRGB2INT(d)
			fits.writeto('tmp3.fits',i[ImgSizeY+NP-1:NP-1:-1,NP:ImgSizeX+NP],overwrite=True)
			#exit()

		edgepixels=lmap!=16777215
		XY=np.argwhere(edgepixels)
		XY0=np.array([self.Edges[lmap[x,y]].p0 for x,y in XY])
		XY1=np.array([self.Edges[lmap[x,y]].p1 for x,y in XY])
		IP0=np.array([self.pmap[x,y] for x,y in XY0])
		IP1=np.array([self.pmap[x,y] for x,y in XY1])
		ok=(XY0[:,0]-XY[:,1])**2+(XY0[:,1]-XY[:,0])**2<(XY1[:,0]-XY[:,1])**2+(XY1[:,1]-XY[:,0])**2
		lmap[XY[ok,0],XY[ok,1]]=IP0[ok]
		lmap[XY[~ok,0],XY[~ok,1]]=IP1[~ok]
		self.filledmap[edgepixels]=lmap[edgepixels]
		if self.pID is not None: self.filledmap=self.pID[self.filledmap-1]
	#END_OF_def CalDVD(self,**kwargs):

	def CalPVD(self,**kwargs):
		for k in self.Edges:
			#print self.Edges[k]
			self.drawLine(self.Edges[k])
		if Voronoi.debug: fits.writeto(self.FileName+'_shell.fits',self.pmap,overwrite=True)
		#Now each cell is enclosed with pixels marking its cell number.
		#Then each cell can be filled using these edge pixels as signposts.
		saveequdist=np.where(self.pmap==-1)
		for x in np.arange(self.ImgSizeX):
			for y in np.arange(self.ImgSizeY):
				if self.pmap[x,y]>0:
					self.pmap[x,0:y] = self.pmap[x,y]
					break
			#My one-direction method fails in case of very large cell which runs through the whole image in the Y-direction
			#Mark such columns of pixels with -2
			if self.pmap[x,0]==0: self.pmap[x,0]=-2
		for x in np.arange(self.ImgSizeX):
			Ncurrent=self.pmap[x,0]
			for y in np.arange(1,self.ImgSizeY):
				if self.pmap[x,y]<=0: self.pmap[x,y]=Ncurrent
				else:
					if self.pmap[x,y]!=Ncurrent: Ncurrent = self.pmap[x,y]
		if Voronoi.debug: fits.writeto(self.FileName+'_face1.fits',self.pmap,overwrite=True)
		if (self.pmap==0).any(): sys.exit('ERROR: sth left')
		self.pmap[saveequdist[0],saveequdist[1]] = -1
		#Now cells are filled, leving only those equidistance pixels which are marked as -1

		#Assign whole column pixels which were marked as -2 to cells
		if np.any(self.pmap==-2):
			xcol,ycol=np.where(self.pmap==-2)
			assert len(np.unique(xcol))*self.ImgSizeY==len(ycol)
			possiblesite=np.unique(self.pmap[xcol-1,:])
			possiblesite=possiblesite[possiblesite>0]
			for x in np.unique(xcol):
				xds=[(x-self.Px[n-1])**2+(y-self.Py[n-1])**2 for n in possiblesite]
				n=possiblesite[np.argmin(xds)]
				self.pmap[x,:]=n
			if Voronoi.debug: fits.writeto(self.FileName+'_face2.fits',self.pmap,overwrite=True)

		#Assign equidistance pixels which were marked as -1 to cells
		VeryRareCase=[] # big hole (some pixels are surrounded by -1 pixels)
		for x,y in np.argwhere(self.pmap==-1):
			Ns = np.array([[fx,fy] for [fx,fy] in [[x+1,y],[x,y+1],[x-1,y],[x,y-1]] \
				if fx>=0 and fy>=0 and fx<self.ImgSizeX and fy<self.ImgSizeY and self.pmap[fx,fy]>0])
			if Ns.size==0:
				Ns = np.array([[fx,fy] for [fx,fy] in [[x+1,y+1],[x+1,y-1],[x-1,y+1],[x-1,y-1]] \
					if fx>=0 and fy>=0 and fx<self.ImgSizeX and fy<self.ImgSizeY and self.pmap[fx,fy]>0])
				if Ns.size==0:
					#very rare case, with very low probability to exist
					VeryRareCase.append([x,y])
					self.pmap[x,y] = 0
					continue
			Ns = self.pmap[Ns[:,0],Ns[:,1]]
			self.pmap[x,y] = 0
			for N in Ns:
				if (Ns==N).sum()>1: self.pmap[x,y] = -N
			if self.pmap[x,y]==0: self.pmap[x,y] = -Ns[np.random.randint(0,Ns.size)]
		self.pmap[self.pmap<0] *= -1
		if Voronoi.debug: fits.writeto(self.FileName+'_face3.fits',self.pmap,overwrite=True)

		for x,y in VeryRareCase:
			Ns = np.array([[fx,fy] for [fx,fy] in [[x+1,y],[x,y+1],[x-1,y],[x,y-1],[x+1,y+1],[x+1,y-1],[x-1,y+1],[x-1,y-1]] \
				if fx>=0 and fy>=0 and fx<self.ImgSizeX and fy<self.ImgSizeY and self.pmap[fx,fy]>0])
			if Ns.size==0: sys.exit('ERROR: pixel %d,%d has no neighbor' %(y+1,x+1))
			Ns = self.pmap[Ns[:,0],Ns[:,1]]
			self.pmap[x,y] = 0
			for N in Ns:
				if (Ns==N).sum()>1: self.pmap[x,y] = N
			if self.pmap[x,y]==0: self.pmap[x,y] = Ns[np.random.randint(0,Ns.size)]
		if (self.pmap==0).any(): sys.exit('ERROR: CalPVD not completed')

	#END_OF_def CalPVD(self)

	def CalAdj(self):
		if Voronoi.debug: print(color("\nCalculating Adjacency List",32,0))
		assert (self.pmap[self.Px,self.Py] == np.arange(1,self.Px.size+1)).all()
		self.Adj = [[] for _ in itertools.repeat(None,self.Px.size+1)]
		for e in self.Edges.values():
			self.Adj[self.pmap[tuple(e.p0)]].append(self.pmap[tuple(e.p1)])
			self.Adj[self.pmap[tuple(e.p1)]].append(self.pmap[tuple(e.p0)])
		#for _ in self.Adj[1:]: assert _
	#END_OF_def CalAdj(self):

