#!/usr/bin/env python
################################################################################
#SweepLine.py
#Voronoi Tessellation with Fortune's sweep line algorithm
#A discrete Voronoi diagram can also be constructed.
#
#Copyright (C) 2014 Teng Liu
#
#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
#
#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#Teng Liu's email: lewtonstein@gmail.com
################################################################################
import numpy as np
import pyfits as pf
import heapq
import itertools
import time,sys,warnings,os,getopt
#warnings.simplefilter('ignore')
color = lambda s,ncolor,nfont: "\033["+str(nfont)+";"+str(ncolor)+"m"+s+"\033[0;0m"
#define:
#p<q if py<qy or py=qy and px<qx
#vertical,upward: y++
#horizontal,right: x++
class HalfEdge(object):
	def __init__(self,p0,p1,direct,base,twin=None,summit=None):
		#if p0[1]<p1[1] or p0[1]==p1[1] and p0[0]<p1[0]:
		#	#the 2 points on 2 sides, p0 < p1
		self.p0 = p0
		self.p1 = p1
		self.direct = direct # 0:upward -1:left 1:right
		self.base = base
		#history test, to save time
		if self.direct == 0:
			self.hist = (self.p0[0]+self.p1[0])/2.
		else:
			self.hist = [None,None]
		self.twin = twin
		self.summit = summit
	def __str__(self):
		return 'p0:'+str(self.p0)+' p1:'+str(self.p1)+' dir:'+str(self.direct)+' b:'+str(self.base)+' s:'+str(self.summit)+' hist:'+str(self.hist)
	def copy(self):
		return HalfEdge(self.p0,self.p1,self.direct,self.base,twin=self.twin,summit=self.summit)
	def isRightOf(self,px,y0):
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
		if self.direct == 0: return self.hist>px
		else:
			if self.hist[0] == y0: cx = self.hist[1]
			else:
				if self.p1[1]>self.p0[1]: u,v = self.p1,self.p0
				else: v,u = self.p1,self.p0
				if u[1]>y0: sys.exit('Error u1>y0   %.8f %.8f - %.8f %.8f - %.8f %.8f' % (u[0],u[1],v[0],v[1],px,y0))
				if round(u[1],Voronoi.SLVround)==round(v[1],Voronoi.SLVround): cx = (u[0]+v[0])/2.
				elif round(u[1],Voronoi.SLVround)==round(y0,Voronoi.SLVround): cx = u[0]
				else:
					insqrt = (u[1]*v[1]-u[1]*y0-v[1]*y0+y0*y0)*(u[0]*u[0]-2*u[0]*v[0]+u[1]*u[1]-2*u[1]*v[1]+v[0]*v[0]+v[1]*v[1])
					cx = (self.direct*np.sqrt(insqrt)-u[0]*v[1]+u[0]*y0+u[1]*v[0]-v[0]*y0)/(u[1]-v[1])
					if insqrt<0: print 'Error sqrt<0 %g %.8f %.8f - %.8f %.8f - %.8f %.8f  %.8f' % (insqrt,u[0],u[1],v[0],v[1],px,y0,cx)
				self.hist = [y0,cx]
			return cx>px
	def Intersect(self,Next):
		if self.p1!=Next.p0: sys.exit('ERROR: not consecutive')
		if self.p0==(-1,-1): return None
		if Next.p1==(-2,-2): return None
		if self.direct<=0 and Next.direct>=0: return None #not Intersect
		#1: a<0 b>0  2: a<0 b=0  3: a=0 b>0  4: a=0 b=0
		ux,uy = self.p0
		vx,vy = self.p1
		wx,wy = Next.p1
		a = (ux-vx)*(vy-wy)-(uy-vy)*(vx-wx)
		if a<=0: return None
		b1 = (ux-vx)*(ux+vx)+(uy-vy)*(uy+vy)
		b2 = (vx-wx)*(vx+wx)+(vy-wy)*(vy+wy)
		cx = (b1*(vy-wy)-b2*(uy-vy))/2./a
		cy = (b2*(ux-vx)-b1*(vx-wx))/2./a
		d = np.sqrt((cx-ux)**2+(cy-uy)**2)
		return cx,cy+d
	def complete(self,xs,ys):
		#fill self.summit
		xb,yb = self.base
		#if xs>xb and not self.direct>0: print self,xs,ys
		#if xs==xb and not self.direct==0: print self,xs,ys
		#if xs<xb and not self.direct<0: print self,xs,ys
		#The directs are good due to the checking here
		if xb<=-0.5 and xs<=-0.5 or xb>=Voronoi.ImgSizeX-0.5 and xs>=Voronoi.ImgSizeX-0.5: self.summit = None
		elif yb<=-0.5 and ys<=-0.5 or yb>=Voronoi.ImgSizeY-0.5 and ys>=Voronoi.ImgSizeY-0.5: self.summit = None
		#elif round(xs,Voronoi.SLVround)==round(xb,Voronoi.SLVround) and round(ys,Voronoi.SLVround)==round(yb,Voronoi.SLVround): self.summit = None #seems useless
		else:
			if not -0.5<=xb<=Voronoi.ImgSizeX-0.5: print 'Unexpected base x',xb #I don't know why x is so good
			if ys>-0.5 and yb<-0.5: #uppward crossing y==0 (bottom)
				self.base = xb-(yb+0.5)*(xs-xb)/(ys-yb), -0.5
				xb,yb = self.base
			elif yb<Voronoi.ImgSizeY-0.5 and ys>Voronoi.ImgSizeY-0.5: #crossing y==max
				#print xb,yb,xs,ys,'T->'
				xs,ys =  xb+(Voronoi.ImgSizeY-0.5-yb)*(xs-xb)/(ys-yb), Voronoi.ImgSizeY-0.5
				#print xb,yb,xs,ys
			elif ys<-0.5 and yb>-0.5: #back crossing y==0
				xs,ys =  xb-(yb+0.5)*(xs-xb)/(ys-yb), -0.5
			elif yb>Voronoi.ImgSizeY-0.5 and ys<Voronoi.ImgSizeY-0.5: #back crossing y==max
				#print xb,yb,xs,ys,'T<-'
				self.base = xb+(Voronoi.ImgSizeY-0.5-yb)*(xs-xb)/(ys-yb), Voronoi.ImgSizeY-0.5
				xb,yb = self.base
				#print xb,yb,xs,ys
			#treat x y in two steps, to cover such cases that crossing two image edges
			if xs<-0.5 and xb>-0.5:
				xs,ys =  -0.5, yb-(xb+0.5)*(ys-yb)/(xs-xb)
			elif xs>Voronoi.ImgSizeX-0.5 and xb<Voronoi.ImgSizeX-0.5:
				#print xb,yb,xs,ys,'R->',
				xs,ys =  Voronoi.ImgSizeX-0.5, yb+(Voronoi.ImgSizeX-0.5-xb)*(ys-yb)/(xs-xb)
				#print xb,yb,xs,ys
			else:
				if not (-0.5<=yb<=Voronoi.ImgSizeY-0.5 and -0.5<=xb<=Voronoi.ImgSizeX-0.5): print 'What',xb,yb,xs,ys
			self.summit = xs,ys
	def line(self,x,y):
		xb,yb = self.base
		if yb<=-0.5 and y<=-0.5 or yb>=Voronoi.ImgSizeY-0.5 and y>=Voronoi.ImgSizeY-0.5:
			return ""
		if y>-0.5 and yb<-0.5:
			#for the edge with (x1,y1) and (x2,y2) on its two sides
			#find the point (x,0) which is equi-distance to the two points
			x1,y1 = self.p0
			x2,y2 = self.p1
			xb = (x2**2 + (y2+0.5)**2 - x1**2 - (y1+0.5)**2) /2. /(x2-x1)
			yb = -0.5
		if yb<Voronoi.ImgSizeY-0.5 and y>Voronoi.ImgSizeY-0.5:
			x1,y1 = self.p0
			x2,y2 = self.p1
			x = (x2**2 + (y2+0.5-Voronoi.ImgSizeY)**2 - x1**2 - (y1+0.5-Voronoi.ImgSizeY)**2) /2. /(x2-x1)
			y = Voronoi.ImgSizeY-0.5
		if y<-0.5 and yb>-0.5:
			x1,y1 = self.p0
			x2,y2 = self.p1
			x = (x2**2 + (y2+0.5)**2 - x1**2 - (y1+0.5)**2) /2. /(x2-x1)
			y = -0.5
		return "line(%.3f,%.3f,%.3f,%.3f) # line=0 0 tag={%d,%d-%d,%d-%d} font=\"helvetica 1 normal roman\"" % (yb+1,xb+1,y+1,x+1,self.p0[0],self.p0[1],self.p1[0],self.p1[1],self.direct)
class SweepTable(dict):
	#a dict of half edges
	def __init__(self,Px,Py):
		self.KL = [] #key list
		self.ysweep=0
		self.p = None

		ymin = Py>Py.min()
		self.Q = EventQueue(Px[ymin],Py[ymin])

		#initialize self dict
		ymin = np.where(Py==Py.min())[0]
		x0,y0 = Px[ymin[0]],Py[ymin[0]]
		self[-1] = HalfEdge((-1,-1),(x0,y0),-1,(-1,-1)) #0: left boundary of image
		self.KL.append(-1)
		self.counter = itertools.count(1) #to make sure each key is different
		for n in np.arange(1,ymin.size):
			x = Px[ymin[n]]
			y = Py[ymin[n]]
			count = self.counter.next()
			self[count] = HalfEdge((x0,y0),(x,y),0,((x0+x)/2.,-0.5))
			self.KL.append(count)
			x0,y0 = x,y
		self[-2] = HalfEdge((x0,y0),(-2,-2),1,(-2,-2)) #1: right boundary of image
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
			Nmid = (Nlow+Nhigh)/2
			#print 'mid:',Nmid
			if self[self.KL[Nmid]].isRightOf(p[0],p[1]):
				Nhigh = Nmid
			else:
				Nlow = Nmid
		return Nlow

	def insert(self,n,x,y):
		#input: self[KL[n]].X <= p < self[KL[n+1]].X
		#p: y,x,...
		x0,y0 = self[self.KL[n]].p1
		base = x,(y**2-y0**2-(x-x0)**2)/2./(y-y0)
		count = self.counter.next()
		self.KL.insert(n+1,count)
		self[count] = HalfEdge((x0,y0),(x,y),-1,base,count+1) #new minus half edge
		#if Voronoi.debug: print 'new: ',count,self[count]
		count = self.counter.next()
		self.KL.insert(n+2,count)
		self[count] = HalfEdge((x,y),(x0,y0),1,base,count-1)  #new plus half edge
		#if Voronoi.debug: print 'new: ',count,self[count]

	def shrink(self):
		x,y,l,r = self.p
		n = self.KL.index(l)
		if r!=self.KL[n+1]: sys.exit('l-not-r')
		else:
			L = self.KL[n-1]
			R = self.KL[n+2]
		self.Q.deleteifhas(L)
		self.Q.deleteifhas(r)
		#(x,y) is the top point of the circle which goes through the three points (l.p0,l.p1,r.p1), l.p1==r.p0
		#the center of the circle is (x,y-r)
		#let x0,y0 be anyone of the three points
		#(x-x0)^2 + (y-r-y0)^2 = r^2
		if self[l].p1!=self[r].p0: sys.exit('no')
		if self[l].p0[1]==self[r].p1[1]:
			direct = 0
		else:
			if self[l].p0[1]>self[r].p1[1]:
				x0 = self[l].p0[0]
				if x==x0: direct = 1
			else:
				x0 = self[r].p1[0]
				if x==x0: direct = -1
			if x < x0: direct = -1
			elif x > x0: direct = 1
		x0,y0 = self[l].p0
		if y==y0: x0,y0 = self[l].p1
		if y==y0: x0,y0 = self[r].p1
		xb,yb = x,(y**2-y0**2-(x-x0)**2)/2./(y-y0)
		#if Voronoi.debug: print 'shrink:',self[l],self[r],x,y,(xb,yb)
		count = self.counter.next()
		self[count] = HalfEdge(self[l].p0,self[r].p1,direct,(xb,yb))
		self.KL[n] = count
		self.KL.pop(n+1)
		self[l].complete(xb,yb)
		self[r].complete(xb,yb)
		i = self[L].Intersect(self[count])
		if i is not None:
			#if Voronoi.debug: print '->V',i,self[L],self[count]
			self.Q.insert(i,L,count)
		i = self[count].Intersect(self[R])
		if i is not None:
			#if Voronoi.debug: print '->V',i,self[count],self[R]
			self.Q.insert(i,count,R)
		return l,r
	def renewEnds(self):
		y = self.p[1]
		todo = []
		if y <= self.ysweep or len(self)<=2: return todo
		if self[self.KL[-2]].isRightOf(Voronoi.ImgSizeX-0.5,y):
			#if Voronoi.debug: print 'EdgeR:',self[-2],self[self.KL[-2]]
			if self[-2].p0!=self[self.KL[-2]].p1: sys.exit('should not')
			u0,u1 = self[self.KL[-2]].p0
			v0,v1 = self[self.KL[-2]].p1
			if u1==v1: sys.exit('cannot')
			yb = (u0**2-2*u0*(Voronoi.ImgSizeX-0.5)+u1**2-v0**2+2*v0*(Voronoi.ImgSizeX-0.5)-v1**2)/2./(u1-v1)
			self[self.KL[-2]].complete(Voronoi.ImgSizeX-0.5,yb)
			todo.append([self.KL[-2],self[self.KL[-2]].copy()])
			if self.KL[-2] in self.Q:
				print 'Please check'
				print self[self.KL[-2]]
				print self[self.KL[-1]]
				exit()
			#maybe the VertexEvt to delete is just the p which has just been deleted from Q as the minimum.
			if self.p is not None and len(self.p)==4 and self.KL[-3]==self.p[2]: self.p = None
			else: self.Q.deleteifhas(self.KL[-3])
			self[-2] = self.pop(self.KL[-2])
			self.KL.pop(-2)
			self[-2].p1 = (-2,-2)
			self[-2].direct = 1
			#if Voronoi.debug: print '-->',self[-2]
		if not self[self.KL[1]].isRightOf(-0.5,y):
			#if Voronoi.debug: print 'EdgeL:',self[-1],self[self.KL[1]]
			if self[-1].p1!=self[self.KL[1]].p0: sys.exit('should not')
			u0,u1 = self[self.KL[1]].p0
			v0,v1 = self[self.KL[1]].p1
			if u1==v1: sys.exit('cannot')
			yb = (u0**2-2*u0*(-0.5)+u1**2-v0**2+2*v0*(-0.5)-v1**2)/2./(u1-v1)
			self[self.KL[1]].complete(-0.5,yb)
			todo.append([self.KL[1],self[self.KL[1]].copy()])
			if self.p is not None and len(self.p)==4 and self.KL[1]==self.p[2]: self.p = None
			else: self.Q.deleteifhas(self.KL[1])
			self[-1] = self.pop(self.KL[1])
			self.KL.pop(1)
			self[-1].p0 = (-1,-1)
			self[-1].direct = -1
			#if Voronoi.debug: print '-->:',self[-1]
		self.ysweep = y
		return todo
class EventQueue(dict):
	def __init__(self,Px,Py):
		if Px.shape==Py.shape==(Px.size,):
			self.S = np.zeros((Px.size,2))
			self.S[:,0] = Py
			self.S[:,1] = Px
			self.SMin = 0
			self.counter = itertools.count(1)
			self.VerEvt = []
		else:
			sys.exit('ERROR: shape error')
	def __str__(self):
		return self.__repr__()+'\n'+str(self.VerEvt)
	def delMin(self):
		while len(self.VerEvt)>0 and self.VerEvt[0][2] == 0: heapq.heappop(self.VerEvt)
		#The interesting property of a heap is that its smallest element is always the root, heap[0].
		if self.SMin<self.S.shape[0] and (len(self.VerEvt)==0 or self.S[self.SMin].tolist() <= self.VerEvt[0][:2]):
			self.SMin += 1
			return self.S[self.SMin-1][[1,0]]	#Site event
		elif len(self.VerEvt)>0 and (self.SMin>=self.S.shape[0] or self.S[self.SMin].tolist() > self.VerEvt[0][:2]):
			y,x,count,l,r = heapq.heappop(self.VerEvt)
			self.pop(l)
			return [x,y,l,r]	#Vertex event
		else:
			return None
		#END_OF_def delMin(self):
	def insert(self,p,l,r,count=None):
		if l in self: sys.exit(str(l)+'already in')
		if count is None: count=self.counter.next()
		evt = [p[1],p[0],count,l,r]
		self[l] = evt
		heapq.heappush(self.VerEvt,evt)
	def deleteifhas(self,l):
		#if Voronoi.debug: print '><:',l
		evt = self.pop(l,None)
		if evt is not None: evt[2] = 0
class Voronoi(object):
	ImgSizeX = 0
	ImgSizeY = 0
	SLVround = 6
	#if SLVround is low, two points far from each other can be taken as neighbors
	debug = False
	#debug = True
	def __init__(self,image=None,events=None,**kwargs):
		#StartTime=time.time()
		self.FileName = kwargs.pop('FileName','FileName')
		print color('Voronoi Construction: '+self.FileName,34,1)
		ToCalArea = kwargs.pop('calarea',False)
		ToCalPVD = kwargs.pop('calpvd',False)
		ToCalAdj = kwargs.pop('caladj',False)
		ToCalDst = kwargs.get('caldst',False)
		Hdr = kwargs.get('Hdr',None)
		self.OffSetX=0
		self.OffSetY=0
		if image is not None: #integer position
			assert events is None
			Voronoi.ImgSizeX,Voronoi.ImgSizeY = image.shape
			self.pmap = np.zeros(image.shape,dtype='int32') #pvd
			self.Py,self.Px = np.where(image.T>0) #coordinates (y,x)
			self.pmap[self.Px,self.Py] = np.arange(1,self.Px.size+1)
		if events is not None: #float position
			assert image is None and events.shape[1] >= 2
			if ToCalPVD: sys.exit('--calpvd not supported')
			if ToCalDst: sys.exit('--caldst not supported')
			if np.min(events[:,0])<0:
			#To be consistent. in integer position case (0~size-1), I used -0.5,size-0.5 as the edge, which corresponds to the ds9 image edge of 0.5 ~ size+0.5
				self.OffSetX = -np.min(events[:,0])
				events[:,0] += self.OffSetX
			if np.min(events[:,1])<0:
				self.OffSetY = -np.min(events[:,1])
				events[:,1] += self.OffSetY
			Voronoi.ImgSizeX,Voronoi.ImgSizeY = kwargs.pop('size',(int(np.max(events[:,:2]))+2,int(np.max(events[:,:2]))+2))
			if np.max(events[:,0])>Voronoi.ImgSizeX-1 or np.max(events[:,1])>Voronoi.ImgSizeY-1:
				print color("ERROR: points out of (0-%g,0-%g)" %(Voronoi.ImgSizeX,Voronoi.ImgSizeY),31,1)
				exit()
			MakeIntImage=kwargs.pop('MakeIntImage',False)
			if MakeIntImage:
				pmap = np.zeros((Voronoi.ImgSizeX,Voronoi.ImgSizeY),dtype='float64')
				reg=file(self.FileName+'_points.reg','w')
				for x,y in events[:,:2]:
					pmap[int(x+0.5),int(y+0.5)] += 1
					print >>reg,"image;circle(%.2f,%.2f,0.3) # color=red;line(%.2f,%.2f,%d,%d) # color=red" % (y+1,x+1, y+1,x+1, int(y+0.5)+1,int(x+0.5)+1)
				if Hdr is None: Hdr=pf.Header()
				Hdr.update('OffSetX',self.OffSetX)
				Hdr.update('OffSetY',self.OffSetY)
				pf. writeto(self.FileName+'_image.fits',pmap,Hdr,clobber=True)
				print '>> '+self.FileName+'_image.fits'
				print '>> '+self.FileName+'_points.reg'
				reg.close()
			Resolution=kwargs.pop('Resolution',3)
			tmp = {(round(events[n,0],Resolution),round(events[n,1],Resolution)):events[n] for n in range(events.shape[0])}
			if len(tmp)!=events.shape[0]:
				print color("Warning: %d duplicated points deleted" % (events.shape[0]-len(tmp)),31,1)
				events = np.array(tmp.values())
			tmp = np.array(sorted(events,key=lambda x:(x[1],x[0]))).round(Voronoi.SLVround)
			self.Px,self.Py = tmp[:,0],tmp[:,1]
			if events.shape[1] == 3: value = tmp[:,2] #to be continued
			self.pmap={(self.Px[n],self.Py[n]):n+1 for n in np.arange(self.Px.size)}
		self.Edges = {}
		self.Amap = None
		self.PPdis = None
		self.Adj = None
		self.EdgePoint = None
		self.dmap = None
		if ToCalDst:
			self.dmap = np.zeros((Voronoi.ImgSizeX,Voronoi.ImgSizeY),dtype='float64')
			self.dmap[:,:] = image[:,:]
			ToCalArea = True # dependent on Amap
			ToCalPVD = True # corrected while CalPVD

		self.Construct(**kwargs)
		if ToCalArea: self.CalArea(**kwargs)
		if ToCalPVD: self.CalPVD(**kwargs)
		if ToCalAdj: self.CalAdj()
		if __name__ == '__main__':
			if ToCalDst:
				print '>> '+self.FileName+'_dst.fits'
				pf.writeto(self.FileName+'_dst.fits',self.dmap,Hdr,clobber=True)
				return
			if ToCalPVD:
				print '>> '+self.FileName+'_pvd.fits'
				pf.writeto(self.FileName+'_pvd.fits',self.pmap,Hdr,clobber=True)
			if ToCalArea:
				if type(self.Amap) is dict:
					print '>> '+self.FileName+'_area.dat'
					#np.savetxt(self.FileName+'_area.dat',self.Amap.values())
					np.savetxt(self.FileName+'_area.dat',np.hstack((np.array(self.Amap.values()).reshape(len(self.Amap),1),np.array(self.Amap.keys())-[self.OffSetX,self.OffSetY])),fmt='%f	%f	%f')
				else:
					print '>> '+self.FileName+'_area.fits'
					pf.writeto(self.FileName+'_area.fits',self.Amap,Hdr,clobber=True)
	#END_OF_init

	def Construct(self,**kwargs):
		StartTime=time.time()
		T = SweepTable(self.Px,self.Py)
		Q = T.Q
		T.p = Q.delMin()
		while T.p is not None:
			for k,e in T.renewEnds():
				self.Edges[k] = e
			if T.p is None: #p is deleted during renewing ending half edges.
				pass
			elif len(T.p)==2: #Site Event
				#if Voronoi.debug: print 'SiteEvent:',T.p
				n = T.locate(T.p)
				L = T.KL[n]
				R = T.KL[n+1]
				T.insert(n,T.p[0],T.p[1])
				l = T.KL[n+1]
				r = T.KL[n+2]
				if R!= T.KL[n+3]: sys.exit('no')
				Q.deleteifhas(L)
				i1 = T[L].Intersect(T[l])
				i2 = T[r].Intersect(T[R])
				if i1 is not None:
					#if Voronoi.debug: print '->V',i1,T[L],T[l]
					Q.insert(i1,L,l)
				if i2 is not None:
					#if Voronoi.debug: print '->V',i2,T[r],T[R]
					Q.insert(i2,r,R)

			elif len(T.p)==4: #Vertex Event
				if T.p[2] not in T or T.p[3] not in T:
					print 'Useless VertexEvent:',T.p
				else:
					#if Voronoi.debug: print 'VertexEvent:',T.p[0],T.p[1]
					l,r = T.shrink()
					self.Edges[l] = T.pop(l)
					self.Edges[r] = T.pop(r)
			else: sys.exit('error2')
			#if Voronoi.debug: print 'T',T
			#if Voronoi.debug: print 'Q',Q
			T.p = Q.delMin()
		#END_OF_while p is not None:
		T.pop(-1)
		T.pop(-2)
		for l in T.keys():
			e = T[l]
			#e.complete((e.p0[0]**2+e.p0[1]**2-e.p1[0]**2-e.p1[1]**2+2*(e.p1[1]-e.p0[1])*(Voronoi.ImgSizeY-0.5))/2/(e.p0[0]-e.p1[0]), Voronoi.ImgSizeY-0.5)
			if e.direct==0: e.complete(e.base[0],Voronoi.ImgSizeY-0.5)
			elif e.direct==1: e.complete(Voronoi.ImgSizeX-0.5, (e.p0[0]**2+e.p0[1]**2-e.p1[0]**2-e.p1[1]**2+2*(e.p1[0]-e.p0[0])*(Voronoi.ImgSizeX-0.5))/2/(e.p0[1]-e.p1[1]))
			elif e.direct==-1: e.complete(-0.5, (e.p0[0]**2+e.p0[1]**2-e.p1[0]**2-e.p1[1]**2-0.25*(e.p1[0]-e.p0[0]))/2/(e.p0[1]-e.p1[1])) #although mostly they cross Voronoi.ImgSizeY, there is a low possibility not. So I input the crossing point at Left or Right edge of image, let complete to recalculate the summit
			self.Edges[l] = T.pop(l)

		for e in self.Edges.values():
			if e.summit is not None:
				if round(e.base[0],Voronoi.SLVround)==round(e.summit[0],Voronoi.SLVround) and round(e.base[1],Voronoi.SLVround)==round(e.summit[1],Voronoi.SLVround):
					e.summit = None
					continue
				if e.twin is not None:
					e2 = self.Edges[e.twin]
					if e2.summit is not None:
						assert e.base==e2.base
						assert self.Edges[e2.twin]==e
						e.base = e2.summit
						e2.summit = None
						assert e.p0==e2.p1 and e.p1==e2.p0
		#Delete nouse edges
		for k in self.Edges.keys():
			if self.Edges[k].summit is None: self.Edges.pop(k)

		#d = np.array([[e.base[1],e.base[0],e.summit[1],e.summit[0]] for e in self.Edges.values() if e.summit is not None])+1
		d = np.array([[e.base[1],e.base[0],e.summit[1],e.summit[0],e.p0[1],e.p0[0],e.p1[1],e.p1[0]] for e in self.Edges.values() if e.summit is not None])+1
		if __name__ == '__main__':
			file(self.FileName+'.reg','w').write('image\n')
			np.savetxt(self.FileName+'.reg',d,fmt="line(%.3f,%.3f,%.3f,%.3f) # tag={%g,%g,%g,%g}")
			print '>> '+self.FileName+'.reg'
			if kwargs.pop('Delaunay',False):
				d = np.array([[e.p0[1],e.p0[0],e.p1[1],e.p1[0],e.base[1],e.base[0],e.summit[1],e.summit[0]] for e in self.Edges.values() if e.summit is not None])+1
				file(self.FileName+'_Delaunay.reg','w').write('image\n')
				np.savetxt(self.FileName+'_Delaunay.reg',d,fmt="line(%g,%g,%g,%g) # tag={%.3f,%.3f,%.3f,%.3f}")
				print '>> '+self.FileName+'_Delaunay.reg'
		if Voronoi.debug: print 'time:',time.time()-StartTime
	#END_OF_Construct()

	def drawLine(self,e):
		if np.isinf(e.summit[0]): print e
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
		if x2==x1: print 'Warning: Zero-Length Edges not eliminated clearly',e
		a = 1.*(y2-y1)/(x2-x1)
		b = y1-a*x1
		x1 = round(x1,Voronoi.SLVround)
		if x1<0: x1 = 0
		elif x1==int(x1): x1 = int(x1)
		else: x1 = int(x1)+1
		x2=int(x2)+1
		if x2>Voronoi.ImgSizeX: x2 = Voronoi.ImgSizeX
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
			if y2<Voronoi.ImgSizeY and self.pmap[X[n],y2]!=-1:
				N = self.pmap[X[n],y2]
				if N>=0 and N<N2: self.pmap[X[n],y2] = N2
		return

	def CalArea(self,**kwargs):
		RemoveEdgePoint = kwargs.pop('RemoveEdgePoint',False)
		ToCalDst = kwargs.pop('caldst',False)
		#for e in self.Edges.values():
		#	assert e.summit is not None
		if Voronoi.debug: print color("\nCalculating Voronoi Cell Areas",32,0)
		if type(self.pmap) is dict:
			self.Amap={(self.Px[n],self.Py[n]):np.float64(0) for n in np.arange(self.Px.size)}
			P0 = np.array([e.p0 for e in self.Edges.values()])
			P1 = np.array([e.p1 for e in self.Edges.values()])
		else:
			self.Amap = np.zeros((Voronoi.ImgSizeX,Voronoi.ImgSizeY),dtype='float64')
			P0 = np.int32([e.p0 for e in self.Edges.values()])
			P1 = np.int32([e.p1 for e in self.Edges.values()])
		E0 = np.array([e.base   for e in self.Edges.values()]).round(Voronoi.SLVround)
		E1 = np.array([e.summit for e in self.Edges.values()]).round(Voronoi.SLVround)

		DisPPList = lambda P0,P1: np.sqrt((P0[:,0]-P1[:,0])**2+(P0[:,1]-P1[:,1])**2)
		self.PPdis = DisPPList(P0,P1)

		#First calculate "Edge area list"
		#"Edge area": area of the triangle made of each edge and one of its according points
		Earea = DisPPList(E0,E1) * self.PPdis / 4.

		#Then calculate Parea
		#"Point area": area of each point (cell)
		for n in np.arange(P0.shape[0]):
			self.Amap[tuple(P0[n])] += Earea[n]
			self.Amap[tuple(P1[n])] += Earea[n]

		#Modify open cells at the image edge
		#Ignore the special cases of the cells at the four corners
		TriArea = lambda p1,p2,p3: abs(p2[0]*p3[1]+p3[0]*p1[1]+p1[0]*p2[1]-p1[0]*p3[1]-p2[0]*p1[1]-p3[0]*p2[1])/2.
		EE = np.vstack((E0,E1))

		########################################################################
		#L R B T [min,max]
		VyL = [Voronoi.ImgSizeY, -1]
		VxB = [Voronoi.ImgSizeX, -1]
		VyR = [Voronoi.ImgSizeY, -1]
		VxT = [Voronoi.ImgSizeX, -1]
		EpL = {}
		EpB = {}
		EpR = {}
		EpT = {}
		EdgePoint=[]
		################################################################################
		#Edge touching (crossing) two image edges
		n = np.where((EE[:,0]==-0.5)&(EE[:,1]==-0.5))[0]	#bottom left
		if n:
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
			EE[n] = [-1,-1]
		n = np.where((EE[:,0]==-0.5)&(EE[:,1]==Voronoi.ImgSizeY-0.5))[0]	#top left
		if n:
			n = n[0]
			x1,y1 = P0[n-P0.shape[0]]
			x2,y2 = P1[n-P0.shape[0]]
			if x1 < x2:
				assert y1 < y2
				EpL[self.pmap[x1,y1]] = [Voronoi.ImgSizeY-0.5]
				EpT[self.pmap[x2,y2]] = [-0.5]
			else:
				assert y1 > y2
				EpL[self.pmap[x2,y2]] = [Voronoi.ImgSizeY-0.5]
				EpT[self.pmap[x1,y1]] = [-0.5]
			EE[n] = [-1,Voronoi.ImgSizeY]
		n = np.where((EE[:,0]==Voronoi.ImgSizeX-0.5)&(EE[:,1]==-0.5))[0]	#bottom right
		if n:
			n = n[0]
			x1,y1 = P0[n-P0.shape[0]]
			x2,y2 = P1[n-P0.shape[0]]
			if x1 < x2:
				assert y1 < y2
				EpB[self.pmap[x1,y1]] = [Voronoi.ImgSizeX-0.5]
				EpR[self.pmap[x2,y2]] = [-0.5]
			else:
				assert y1 > y2
				EpB[self.pmap[x2,y2]] = [Voronoi.ImgSizeX-0.5]
				EpR[self.pmap[x1,y1]] = [-0.5]
			EE[n] = [Voronoi.ImgSizeX,-1]
		n = np.where((EE[:,0]==Voronoi.ImgSizeX-0.5)&(EE[:,1]==Voronoi.ImgSizeY-0.5))[0]	#top right
		if n:
			n = n[0]
			x1,y1 = P0[n-P0.shape[0]]
			x2,y2 = P1[n-P0.shape[0]]
			if x1 < x2:
				assert y1 > y2
				EpT[self.pmap[x1,y1]] = [Voronoi.ImgSizeX-0.5]
				EpR[self.pmap[x2,y2]] = [Voronoi.ImgSizeY-0.5]
			else:
				assert y1 < y2
				EpT[self.pmap[x2,y2]] = [Voronoi.ImgSizeX-0.5]
				EpR[self.pmap[x1,y1]] = [Voronoi.ImgSizeY-0.5]
			EE[n] = [Voronoi.ImgSizeX,Voronoi.ImgSizeY]
		################################################################################
		#Edges touching (crossing) one image edge
		if Voronoi.debug: EPfile = file('EdgePoint.reg','w')
		for n in np.where(EE[:,0]==-0.5)[0]:	#left EpL
			N = self.pmap[tuple(P0[n-P0.shape[0]])]
			EdgePoint.append(N)
			EpL[N] = EpL.get(N,[])+[EE[n,1]]
			N = self.pmap[tuple(P1[n-P0.shape[0]])]
			EdgePoint.append(N)
			EpL[N] = EpL.get(N,[])+[EE[n,1]]
			VyL[0] = min(VyL[0],EE[n,1])
			VyL[1] = max(VyL[1],EE[n,1])
		for n in np.where(EE[:,1]==-0.5)[0]:	#bottom EpB
			N = self.pmap[tuple(P0[n-P0.shape[0]])]
			EdgePoint.append(N)
			EpB[N] = EpB.get(N,[])+[EE[n,0]]
			N = self.pmap[tuple(P1[n-P0.shape[0]])]
			EdgePoint.append(N)
			EpB[N] = EpB.get(N,[])+[EE[n,0]]
			VxB[0] = min(VxB[0],EE[n,0])
			VxB[1] = max(VxB[1],EE[n,0])
		for n in np.where(EE[:,0]==Voronoi.ImgSizeX-0.5)[0]:	#right EpR
			N = self.pmap[tuple(P0[n-P0.shape[0]])]
			EdgePoint.append(N)
			EpR[N] = EpR.get(N,[])+[EE[n,1]]
			if Voronoi.debug: print >>EPfile,'line(%.2f,%.2f,%.2f,%.2f)' %(self.Py[N-1]+1,self.Px[N-1]+1,EE[n,1]+1,EE[n,0]+1)
			N = self.pmap[tuple(P1[n-P0.shape[0]])]
			EdgePoint.append(N)
			EpR[N] = EpR.get(N,[])+[EE[n,1]]
			if Voronoi.debug: print >>EPfile,'line(%.2f,%.2f,%.2f,%.2f)' %(self.Py[N-1]+1,self.Px[N-1]+1,EE[n,1]+1,EE[n,0]+1)
			VyR[0] = min(VyR[0],EE[n,1])
			VyR[1] = max(VyR[1],EE[n,1])
		for n in np.where(EE[:,1]==Voronoi.ImgSizeY-0.5)[0]:	#top EpT
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
				print >>EPfile,"circle(%d,%d,2) # color=red" % (self.Py[N-1]+1,self.Px[N-1]+1)
			EPfile.close()
		########################################################################
		# For such case where there is no crossing point in one edge of the image
		if len(EpL)==0:
			assert VyL == [Voronoi.ImgSizeY, -1]
			n = np.argmin(self.Px)
			N = self.pmap[self.Px[n],self.Py[n]]
			assert len(EpB[N])==1 and len(EpT[N])==1 and EpB[N][0]==VxB[0] and EpT[N][0]==VxT[0]
			EpL[N] = [-0.5, Voronoi.ImgSizeX-0.5]
			EpB[N].append(-0.5)
			EpT[N].append(-0.5)
		if len(EpB)==0:
			assert VxB == [Voronoi.ImgSizeX, -1]
			n = np.argmin(self.Py)
			N = self.pmap[self.Px[n],self.Py[n]]
			assert len(EpL[N])==1 and len(EpR[N])==1 and EpL[N][0]==VyL[0] and EpR[N][0]==VyR[0]
			EpB[N] = [-0.5, Voronoi.ImgSizeY-0.5]
			EpL[N].append(-0.5)
			EpR[N].append(-0.5)
		if len(EpR)==0:
			assert VyR == [Voronoi.ImgSizeY, -1]
			n = np.argmax(self.Px)
			N = self.pmap[self.Px[n],self.Py[n]]
			assert len(EpB[N])==1 and len(EpT[N])==1 and EpB[N][0]==VxB[1] and EpT[N][0]==VxT[1]
			EpR[N] = [-0.5, Voronoi.ImgSizeX-0.5]
			EpB[N].append(Voronoi.ImgSizeX-0.5)
			EpT[N].append(Voronoi.ImgSizeX-0.5)
		if len(EpT)==0:
			assert VxT == [Voronoi.ImgSizeX, -1]
			n = np.argmax(self.Py)
			N = self.pmap[self.Px[n],self.Py[n]]
			assert len(EpL[N])==1 and len(EpR[N])==1 and EpL[N][0]==VyL[1] and EpR[N][0]==VyR[1]
			EpT[N] = [-0.5, Voronoi.ImgSizeY-0.5]
			EpL[N].append(Voronoi.ImgSizeY-0.5)
			EpR[N].append(Voronoi.ImgSizeY-0.5)
		########################################################################
		for N in EpL.keys():
			x,y = self.Px[N-1],self.Py[N-1]
			if len(EpL[N])==1:
				if EpL[N][0]==VyL[0] and N in EpB.keys():
					#N refers to the left bottom pixel
					assert len(EpB[N])==1 and EpB[N][0]==VxB[0]
					EpL[N].append(-0.5)
					EpB[N].append(-0.5)
				elif EpL[N][0]==VyL[1] and N in EpT.keys():
					#Only one cross point in this edge. N refers to the left top pixel
					assert len(EpT[N])==1 and EpT[N][0]==VxT[0]
					EpL[N].append(Voronoi.ImgSizeX-0.5)
					EpT[N].append(-0.5)
				else:
					print 'EpL',x,y,EpL[N]
					raise RuntimeError('ERROR: Edge point Left')
		for N in EpR.keys():
			x,y = self.Px[N-1],self.Py[N-1]
			if len(EpR[N])==1:
				if EpR[N][0]==VyR[0] and N in EpB.keys():
					assert len(EpB[N])==1 and EpB[N][0]==VxB[1]
					EpR[N].append(-0.5)
					EpB[N].append(Voronoi.ImgSizeX-0.5)
				elif EpR[N][0]==VyR[1] and N in EpT.keys():
					assert len(EpT[N])==1 and EpT[N][0]==VxT[1]
					EpR[N].append(Voronoi.ImgSizeX-0.5)
					EpT[N].append(Voronoi.ImgSizeX-0.5)
				else:
					print 'EpR',x,y,EpR[N]
					n= np.where((P0[:,0]==x)&(P0[:,1]==y))
					if len(n)>0:
						n=n[0][0]
						print 'P0',P0[n,0],EE[n],EE[n-P0.shape[0]]
					n= np.where((P1[:,0]==x)&(P1[:,1]==y))
					if len(n)>0:
						n=n[0][0]
						print 'P1',P1[n,0],EE[n],EE[n-P1.shape[0]]
					raise RuntimeError('ERROR: Edge point Right')
		########################################################################
		for N in EpL.keys():
			x,y = self.Px[N-1],self.Py[N-1]
			if len(EpL[N])!=2: print y+1,x+1,EpL[N]
			assert len(EpL[N])==2
			self.Amap[x,y] += (x+0.5)*abs(EpL[N][0]-EpL[N][1])/2.
			#print "circle(%.1f,%.1f,1.3) # text={%.3f}" %(y+1,x+1,(x+0.5)*abs(EpL[N][0]-EpL[N][1])/2.)
		for N in EpB.keys():
			x,y = self.Px[N-1],self.Py[N-1]
			if len(EpB[N])!=2: print y+1,x+1,EpB[N]
			assert len(EpB[N])==2
			self.Amap[x,y] += (y+0.5)*abs(EpB[N][0]-EpB[N][1])/2.
			#print "circle(%.1f,%.1f,1.3) # text={%.3f}" %(y+1,x+1,(y+0.5)*abs(EpB[N][0]-EpB[N][1])/2.)
		for N in EpR.keys():
			x,y = self.Px[N-1],self.Py[N-1]
			if len(EpR[N])!=2: print y+1,x+1,EpR[N]
			assert len(EpR[N])==2
			self.Amap[x,y] += (Voronoi.ImgSizeX-0.5-x)*abs(EpR[N][0]-EpR[N][1])/2.
			#print "circle(%.1f,%.1f,1.3) # text={%.3f}" %(y+1,x+1,(Voronoi.ImgSizeX-0.5-x)*abs(EpR[N][0]-EpR[N][1])/2.)
		for N in EpT.keys():
			x,y = self.Px[N-1],self.Py[N-1]
			if len(EpT[N])!=2: print y+1,x+1,EpT[N],N
			assert len(EpT[N])==2
			self.Amap[x,y] += (Voronoi.ImgSizeY-0.5-y)*abs(EpT[N][0]-EpT[N][1])/2.
			#print "circle(%.1f,%.1f,1.3) # text={%.3f}" %(y+1,x+1,(Voronoi.ImgSizeY-0.5-y)*abs(EpT[N][0]-EpT[N][1])/2.)
		########################################################################
		if ToCalDst:
			assert np.all(self.Amap[self.Px,self.Py]>0)
			self.dmap[self.Px,self.Py] /= self.Amap[self.Px,self.Py]
		if RemoveEdgePoint:
			for N in self.EdgePoint:
				self.Amap[self.Px[N-1],self.Py[N-1]] = 0
			print color('EdgePoints removed in Amap',34,1)
	#END_OF_def CalArea(self):

	def CalPVD(self,**kwargs):
		ToCalDst = kwargs.pop('caldst',False)

		for k in self.Edges.keys():
			#print self.Edges[k]
			self.drawLine(self.Edges[k])
		if Voronoi.debug: pf.writeto(self.FileName+'_shell.fits',self.pmap,clobber=True)
		#Now each cell is enclosed with pixels marking its cell number.
		#Then each cell can be filled using these edge pixels as signposts.

		##I don't underestand these old codes. Anyway the new codes below works well.
		##Fill self.pmap using edge lines as marks
		#for x in np.arange(Voronoi.ImgSizeX):
		#	print x,'------------------------'
		#	print self.pmap[x,:]
		#	if self.pmap[x,0]==0 or self.pmap[x,1]==0:
		#		y0 = 0
		#		N = self.pmap[x,0]
		#		setNew = False
		#	else:
		#		setNew = True
		#	for y in np.arange(1,Voronoi.ImgSizeY):
		#		if self.pmap[x,y] == -1:
		#			if not setNew:
		#				if N==0: print "Can't do anything"
		#				self.pmap[x,y0:y] = N
		#				setNew = True
		#		elif self.pmap[x,y] > 0:
		#			print y,self.pmap[x,y],setNew
		#			if not setNew:
		#				if y+1<Voronoi.ImgSizeY and self.pmap[x,y+1]!=self.pmap[x,y] and self.pmap[x,y+1]!=0:
		#					if N==0: N = self.pmap[x,y]
		#					self.pmap[x,y0:y] = N
		#					setNew = True
		#			else:
		#				if y+1<Voronoi.ImgSizeY and self.pmap[x,y+1]==0:
		#					y0 = y
		#					N = self.pmap[x,y]
		#					setNew = False
		#	if not setNew:
		#		if Voronoi.debug: print x,y0,'~',y,'->',N
		#		self.pmap[x,y0:y+1] = N
		saveequdist=np.where(self.pmap==-1)
		for x in np.arange(Voronoi.ImgSizeX):
			for y in np.arange(Voronoi.ImgSizeY):
				if self.pmap[x,y]>0:
					self.pmap[x,0:y] = self.pmap[x,y]
					break
			assert self.pmap[x,0]>0
		for x in np.arange(Voronoi.ImgSizeX):
			Ncurrent=self.pmap[x,0]
			for y in np.arange(1,Voronoi.ImgSizeY):
				if self.pmap[x,y]<=0: self.pmap[x,y]=Ncurrent
				else:
					if self.pmap[x,y]!=Ncurrent: Ncurrent = self.pmap[x,y]
		if Voronoi.debug: pf.writeto(self.FileName+'_line1.fits',self.pmap,clobber=True)
		if (self.pmap==0).any(): sys.exit('ERROR: sth left')
		self.pmap[saveequdist[0],saveequdist[1]] = -1
		#Now cells are filled, leving only those equidistance pixels which are marked as -1

		if ToCalDst:
			self.dmap[self.pmap>0] = self.dmap[self.Px,self.Py][self.pmap[self.pmap>0]-1]
			#vmap is not uniquely defined, thus by now the uncertainty is transfered into dmap

		#Assign equidistance pixels which are marked as -1 to cells
		VeryRareCase=[] # big hole
		for x,y in np.argwhere(self.pmap==-1):
			Ns = np.array([[fx,fy] for [fx,fy] in [[x+1,y],[x,y+1],[x-1,y],[x,y-1]] \
				if fx>=0 and fy>=0 and fx<Voronoi.ImgSizeX and fy<Voronoi.ImgSizeY and self.pmap[fx,fy]>0])
			if Ns.size==0:
				Ns = np.array([[fx,fy] for [fx,fy] in [[x+1,y+1],[x+1,y-1],[x-1,y+1],[x-1,y-1]] \
					if fx>=0 and fy>=0 and fx<Voronoi.ImgSizeX and fy<Voronoi.ImgSizeY and self.pmap[fx,fy]>0])
				if Ns.size==0:
					#very rare case, with very low probability to exist
					VeryRareCase.append([x,y])
					self.pmap[x,y] = 0
					continue
			if ToCalDst: self.dmap[x,y] = np.mean(self.dmap[Ns[:,0],Ns[:,1]])
			Ns = self.pmap[Ns[:,0],Ns[:,1]]
			self.pmap[x,y] = 0
			for N in Ns:
				if (Ns==N).sum()>1: self.pmap[x,y] = -N
			if self.pmap[x,y]==0: self.pmap[x,y] = -Ns[np.random.randint(0,Ns.size)]
		self.pmap[self.pmap<0] *= -1
		#Done

		for x,y in VeryRareCase:
			Ns = np.array([[fx,fy] for [fx,fy] in [[x+1,y],[x,y+1],[x-1,y],[x,y-1],[x+1,y+1],[x+1,y-1],[x-1,y+1],[x-1,y-1]] \
				if fx>=0 and fy>=0 and fx<Voronoi.ImgSizeX and fy<Voronoi.ImgSizeY and self.pmap[fx,fy]>0])
			if Ns.size==0: sys.exit('ERROR: pixel %d,%d has no neighbor' %(y+1,x+1))
			if ToCalDst: self.dmap[x,y] = np.mean(self.dmap[Ns[:,0],Ns[:,1]])
			Ns = self.pmap[Ns[:,0],Ns[:,1]]
			self.pmap[x,y] = 0
			for N in Ns:
				if (Ns==N).sum()>1: self.pmap[x,y] = N
			if self.pmap[x,y]==0: self.pmap[x,y] = Ns[np.random.randint(0,Ns.size)]
		if (self.pmap==0).any(): sys.exit('ERROR: CalPVD not completed')
	#END_OF_def CalPVD(self)

	def CalAdj(self):
		if Voronoi.debug: print color("\nCalculating Adjacency List",32,0)
		assert (self.pmap[self.Px,self.Py] == np.arange(1,self.Px.size+1)).all()
		self.Adj = [[] for _ in itertools.repeat(None,self.Px.size+1)]
		for e in self.Edges.values():
			self.Adj[self.pmap[tuple(e.p0)]].append(self.pmap[tuple(e.p1)])
			self.Adj[self.pmap[tuple(e.p1)]].append(self.pmap[tuple(e.p0)])
		#for _ in self.Adj[1:]: assert _
	#END_OF_def CalAdj(self):

#1.Useless check is not necessary anylonger, it's useless
#2.locate maybe improved to be faster
##By profile.run(), I found that __str__() for large data structure is very slow, sometimes printing it can be the most time-consuming function.
def main():
	def usage():
		print """SweepLine.py File [options]

The input File can be image (.fits file) or events (ASCII file storing the point coordinates in its first two columns).

OPTIONS
	-D				Make Delaunay diagram.
	--calarea		Calculate cell Areas.
	--calpvd		Calculate Discrete Voronoi Diagram, only in case of image input.
	--caldst		Calculate Density image, only in case of image input.
	--size M,N		Set image size, only in case of events input.
	--resolution N		Set position resolution to the Nth decimial place.
				Only in case of events input, where high resolution can make problems.
	--makeimage		Make an image from events list
	--rmedgepoint		Remove edge points in Area map
	-h/--help		Help

Caveat
	PVD (Discrete Voronoi Diagram) is not uniquely defined. There is uncertainty in some pixels.
		"""
		exit()
	Options={}
	S_opt='dDh'
	L_opt=['calpvd','calarea','caldst','rmedgepoint','size=','resolution=','accuracy=','makeimage','help']
	opts,args=getopt.getopt(sys.argv[1:],S_opt,L_opt)
	if len(args)>0:
		for arg in args:
			if os.path.isfile(arg):
				InputFile=arg
				sys.argv.remove(arg)
				break
		opts,args=getopt.getopt(sys.argv[1:],S_opt,L_opt)
	for opt,arg in opts:
		if opt == '--calpvd':
			Options['calpvd']=True
		elif opt == '-D':
			Options['Delaunay']=True
		elif opt == '--calarea':
			Options['calarea']=True
		elif opt == '--caldst':
			Options['caldst']=True
		elif opt == '--makeimage':
			Options['MakeIntImage']=True
		elif opt == '--rmedgepoint':
			Options['RemoveEdgePoint']=True
		elif opt == '--size':
			try:
				xy = arg.replace('x',',').split(',')
				x,y = int(xy[0]),int(xy[1])
			except:
				print 'ERROR: --size',arg
				exit()
			else:
				Options['size'] = x,y
		elif opt == '--resolution':
			try:
				n = int(arg)
				assert n>=0 and n<Voronoi.SLVround
			except:
				print 'ERROR: --resolution',arg
				exit()
			else:
				Options['Resolution'] = n
		elif opt == '--accuracy':
			try:
				n = int(arg)
				assert n>0
			except:
				print 'ERROR: --accuracy',arg
				exit()
			else:
				Voronoi.SLVround = n
		elif opt == "-h" or opt == "--help":
			usage()
		elif opt == '-d':
			Voronoi.debug = True
	if len(args)>0: sys.exit("I don't understand"+str(args))
	if 'InputFile' not in locals():
		print "Please specify one image!\n"
		usage()

	if len(InputFile.rsplit('.',1))>1 and InputFile.rsplit('.',1)[1] == 'fits':
		Voronoi(image=pf.getdata(InputFile),FileName=InputFile.rsplit('.',1)[0],Hdr=pf.getheader(InputFile),**Options)
	else: #a file which store the coordinates of points in the first two columns
		Voronoi(events=np.loadtxt(InputFile),FileName=InputFile.rsplit('.',1)[0],**Options)
#import profile
if __name__ == '__main__':
	#profile.run('main()')
	main()
