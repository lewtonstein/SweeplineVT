import setuptools
from lib.CweeplineVT import cc

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="SweeplineVT",
    version="1.0.0",
    author="Teng Liu",
    author_email="lewtonstein@gmail.com",
    description="Voronoi Tessellation using Sweep-line algorithm",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/lewtonstein/SweeplineVT",
    packages=setuptools.find_packages(),
	scripts=['bin/pl_VT.py','bin/slvt.py'],
    classifiers=[
        "Programming Language :: Python :: 3",
		"License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
    ],
	install_requires=["numpy", "astropy", "matplotlib", "numba"],
    python_requires='>=3.6',
	ext_modules=[cc.distutils_extension()]
)
