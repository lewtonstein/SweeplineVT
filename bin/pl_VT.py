#!/usr/bin/env python
import os,sys
import numpy as np
import pylab as pl
import json
__doc__="""
Input the *_VT.dat file from slvt.py and optionally:
	the *_ctd.dat file from slvt.py,
	or "-s", which means plotting step by step.
"""

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
	pl.plot([l[1],l[3]],[l[2],l[4]],'b-',lw=1)
	pl.plot([l[5],l[7]],[l[6],l[8]],'ro')
if d.shape[1]==11:
	dist2=d[:,9]**2+d[:,10]**2
	l=d[np.argsort(dist2)[-1]]
	pl.plot([l[1],l[3]],[l[2],l[4]],'-',lw=3)
	pl.plot([l[5],l[7]],[l[6],l[8]],'C2v',ms=8,lw=3)
	l=d[np.argsort(dist2)[-2]]
	pl.plot([l[1],l[3]],[l[2],l[4]],'-',lw=3)
	pl.plot([l[5],l[7]],[l[6],l[8]],'C2+',ms=8,lw=3)
	l=d[np.argsort(dist2)[-3]]
	pl.plot([l[1],l[3]],[l[2],l[4]],'-',lw=3)
	pl.plot([l[5],l[7]],[l[6],l[8]],'C2x',ms=8,lw=3)
if Step: pl.ioff()
if len(sys.argv)>2 and os.path.isfile(sys.argv[2]): 
	ctd=np.loadtxt(sys.argv[2])
	for l in ctd:
		pl.plot([l[1],l[3]],[l[2],l[4]],'g-')
pl.savefig(VTfile.rsplit('.')[0]+'.png',bbox_inches='tight')
print('>> '+VTfile.rsplit('.')[0]+'.png')
pl.tight_layout()
pl.show()
