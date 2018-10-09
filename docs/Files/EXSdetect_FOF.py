#!/usr/bin/env python
import numpy as np
import sys,itertools
def Percolate(Nplus=1,pmap=None,Connections=None,Vertices=None,root=None,GetCut=False,GetOne=False):
	'''
	Percolate through all the negative pixels in pmap
	Grouping them into sources, each of which is marked with an unique positive number (Nplus).
	"Nplus" starts from 1, so that it means the number of detected sources +1, or the number of next source.
	"Vertices" contains the indices of the involved points.
	"Connections" are the edges of the graph, saved in the format that Connections[N] is a list of point indices which are connected to the point N.
	'''
	def DFTFlood(pmap,p0,Nplus):
		#pmap:  <0: not visited; ==0: useless; >0: visited
		stack = [p0]
		while stack:
			p = stack.pop()
			if pmap[p]>=0: continue
			pmap[p] = Nplus
			for p in [(p[0]+1,p[1]),(p[0]-1,p[1]),(p[0],p[1]+1),(p[0],p[1]-1)]:
				try:
					if pmap[p]<0: stack.append(p)
				except IndexError,msg:
					pass
				except:
					sys.exit('Unexpected Error')
	def DFTraversal_G(Connections, root, visited, Nplus):
		#Connections is the adjacency list/dict of a graph
		#root is the starting site
		#visited is the marker dict, in which ==0: not visited; >0: visited
		#also visited is used as the input site list, sites not included in it are ignored
		#Nplus is a >0 number
		stack = [(root,root)]
		DFN = itertools.count(1)
		while stack:
			parent, site = stack.pop()
			if visited.get(site,-1)==0:
				visited[site] = Nplus
				yield parent,site,DFN.next()
				stack.extend([(site,child) for child in Connections[site] if visited.get(child,-1)==0])
		#this function is a generator as long as there is a "yield", no matter it is run or not.
	def DFTraversal_noG(Connections, root, visited, Nplus):
		stack = [root]
		while stack:
			site = stack.pop()
			if visited.get(site,-1)==0:
				visited[site] = Nplus
				stack.extend([child for child in Connections[site] if visited.get(child,-1)==0])
	################################################################################	
	if pmap is None:
		assert Connections is not None and Vertices is not None
		visited = dict.fromkeys(Vertices,0) #==0 means not visited
		if not GetCut:
			if GetOne:
				assert root is not None
				DFTraversal_noG(Connections, root, visited, Nplus)
				return visited #GetOne, not GetCut, Connections mode
			if root is not None:
				DFTraversal_noG(Connections, root, visited, Nplus)
				Nplus += 1
			for n in (n for n in Vertices if visited[n]==0):
				DFTraversal_noG(Connections, n, visited, Nplus)
				Nplus += 1
			return Nplus,visited #default, not GetOne, not GetCut, Connections mode
		else: #GetCut
			############################################################
			#starting point is given, search only once and return the search tree
			assert root is not None
			tree = {} #value is the child of key
			N_DFT = {}
			N_low = {}
			queue = []
			cutpoint = []
			for parent,n,DFN in DFTraversal_G(Connections, root, visited, Nplus):
				N_DFT[n] = DFN
				queue.append(n)
				if tree.has_key(parent): tree[parent].append(n)
				else: tree[parent] = [n]
				if not tree.has_key(n): tree[n] = []
			for n in queue[:0:-1]:
				Neighbor = [N_DFT[i] for i in Connections[n] if visited.get(i,0)>0 and n not in tree[i]]
				if Neighbor: N_low[n] = min(Neighbor)
				else: N_low[n] = N_DFT[n]+1
				if tree[n]: #not empty, has child
					if N_DFT[n]<=max([N_low[i] for i in tree[n]]): cutpoint.append(n)
					N_low[n] = min(N_low[n],min([N_low[i] for i in tree[n]]))
			tree[root].remove(root)
			if len(tree[root])>1: cutpoint.append(root)
			return visited,cutpoint
	else: #pmap mode
		#"friend" is defined as the neighboring pixel (distance==1).
		#Here the adjacency list of the graph is not used, the pmap is enough
		for p in (p for p in np.argwhere(pmap<0) if pmap[tuple(p)]<0):
			DFTFlood(pmap, tuple(p), Nplus)
			Nplus += 1
		assert (pmap>=0).all()
		return Nplus
#END_OF_def Percolate()

if __name__=='__main__':
	'''
	Input a file whose first two columns contain the indices of each pair of points which are connected in the sense of, for example, being within a distance.
	Such a file is generally the result of a cross-correlation.
	'''
	try:
		edges=np.loadtxt(sys.argv[1],usecols=(0,1),dtype='int32')
	except:
		sys.exit('ERROR: reading the input file')
	else:
		Connections = {}
		for n1,n2 in edges:
			Connections[n1] = Connections.get(n1,[])+[n2]
			Connections[n2] = Connections.get(n2,[])+[n1]
		V = Connections.keys()
		Nplus,visited = Percolate(Connections=Connections,Vertices=V)
		for N in V:
			print visited[N],N
