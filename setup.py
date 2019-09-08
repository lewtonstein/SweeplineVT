import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="SweeplineVT",
    version="0.0.1",
    author="Teng Liu",
    author_email="lewtonstein@gmail.com",
    description="Voronoi Tessellation using Sweep-line algorithm",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/lewtonstein/SweepLine",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
		"License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
