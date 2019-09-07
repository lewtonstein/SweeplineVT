#!/usr/bin/env python
import os,sys
import numpy as np
import pylab as pl
import json
if '-s' in sys.argv:
	sys.argv.remove('-s')
	Step=True
	pl.ion()
else: Step=False

VTfile=sys.argv[1]
with open(VTfile,'r') as fin:
	l=fin.readline()
	h=json.loads(l.strip('#'))
	d=np.loadtxt(fin)
fig,ax=pl.subplots()
pl.xlim(h['xlow'],h['xhigh'])
pl.ylim(h['ylow'],h['yhigh'])
ax.set(adjustable='box', aspect='equal')
for l in d:
	if Step:
		input('Enter to continue')
		pl.draw()
	pl.plot([l[1],l[3]],[l[2],l[4]],'b-')
	pl.plot([l[5],l[7]],[l[6],l[8]],'ro')
if Step: pl.ioff()
if len(sys.argv)>2 and os.path.isfile(sys.argv[2]): 
	ctd=np.loadtxt(sys.argv[2])
	for l in ctd:
		pl.plot([l[1],l[3]],[l[2],l[4]],'g-')
pl.savefig(VTfile.rsplit('.')[0]+'.png',bbox_inches='tight')
print('>> '+VTfile.rsplit('.')[0]+'.png')
pl.tight_layout()
pl.show()
