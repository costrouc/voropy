"""
Voropy
******

A 3D cell-based Voronoi library based on voro++
-----------------------------------------------

`Code available`_ on Github.

`Full documentation available`_ at Read the Docs.

.. _Code available: https://github.com/costrouc/voropy

.. _Full documentation available: https://voropy.readthedocs.org

Description
-----------

Voropy is a library to calculate Voronoi (and Laguerre) tessellations in 3D and analyze their
structure. The tessellation is calculated as a `list` of `Cell` objects, each of which
can give information on its volume, centroid, number of faces, surface area, etc. The library is
made with packings of spherical particles in mind, possibly with variable sizes.

"""

from setuptools import setup, find_packages
from distutils.extension import Extension
from distutils.command.sdist import sdist as _sdist

try:
    from Cython.Build import cythonize
except ImportError:
    cythonize = None


if cythonize is not None:
    print("Building with Cython.")
    ext = cythonize(Extension("voropy._voro",
              sources=["voropy/_voro.pyx", "src/voro++.cc"],
              include_dirs=["src"],
              language="c++",
    ))
else:
    print("Cython not found, using prebuilt file.")
    ext = [Extension("voropy._voro",
              sources=["voropy/_voro.cpp", "src/voro++.cc"],
              include_dirs=["src"],
              language="c++",
    )]


# Set sdist to make the .cpp file
# from http://stackoverflow.com/a/18418524/4190270
class sdist(_sdist):
    def run(self):
        # this is already imported, but the import might have failed. If so, raise an ImportError now.
        from Cython.Build import cythonize

        # Make sure the compiled Cython files in the distribution are up-to-date
        cythonize(['voropy/_voro.pyx'])
        _sdist.run(self)

cmdclass = dict(sdist = sdist)

# create the extension and add it to the python distribution
setup(
    name='voropy',
    version='0.2',
    author="Chris Ostrouchov",
    author_email="chris.ostrouchov@gmail.com",
    description = ("A module for calculating and analyzing 3D Voronoi tessellations"),
    license = "BSD",
    keywords = "laguerre voronoi tessellation voro++",
    url = "https://voropy.readthedocs.org",
    long_description=__doc__.strip(),
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Topic :: Scientific/Engineering :: Mathematics",
        "Topic :: Scientific/Engineering :: Physics",
        "Intended Audience :: Science/Research",
        "Operating System :: POSIX :: Linux",
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: BSD License",
    ],
    tests_require=['pytest'],
    packages=find_packages(exclude=['tests']),
    ext_modules = ext,
    cmdclass=cmdclass,
)
