import setuptools
import numpy as np
from numba.pycc import CC
cc = CC('CweeplineVT')
cc.verbose = True
@cc.export('CalCxRightOf','float64(float64,float64,float64,float64,int32,float64,float64,float64)')
def CalCxRightOf(u0,u1,v0,v1,self_direct,px,y0,Voronoi_atol):
	if abs(u1-v1)<Voronoi_atol: cx = (u0+v0)/2.  #np.isclose is very slow
	elif abs(u1-y0)<Voronoi_atol: cx = u0
	elif u0==v0==px: cx=np.float64(np.inf)*self_direct #fastmath=True doesn't work with inf
	else:
		insqrt = (u1*v1-u1*y0-v1*y0+y0*y0)*(u0*u0-2*u0*v0+u1*u1-2*u1*v1+v0*v0+v1*v1)
		cx = (self_direct*np.sqrt(insqrt)-u0*v1+u0*y0+u1*v0-v0*y0)/(u1-v1)
		if insqrt<0: print('Error sqrt<0',insqrt,u0,u1,v0,v1,px,y0,cx)
	return cx #has to be sent back to the Python object to set the hist=(y0,cx)
@cc.export('CalIntersect','float64[:](float64,float64,float64,float64,float64,float64)')
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

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="SweeplineVTtest-lewton",
    version="1.0.0",
    author="Teng Liu",
    author_email="lewtonstein@gmail.com",
    description="Voronoi Tessellation using Sweep-line algorithm",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/lewtonstein/SweeplineVT",
	packages=["SweeplineVT"],
	package_dir={"SweeplineVT": "SweeplineVT"},
	scripts=['bin/pl_VT.py','bin/slvt.py'],
    classifiers=[
        "Programming Language :: Python :: 3",
		"License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
    ],
	install_requires=["numpy", "astropy", "matplotlib", "numba","tqdm"],
    python_requires='>=3.6',
	ext_modules=[cc.distutils_extension()]
)
	#package_dir={"": "SweeplineVT"},
    #packages=setuptools.find_packages(where="SweeplineVT"),
    #packages=setuptools.find_packages(),
