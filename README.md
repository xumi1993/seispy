# seispy

[![pipeline status](https://img.shields.io/travis/com/xumi1993/seispy)](https://travis-ci.com/xumi1993/seispy)
[![Build Status](https://img.shields.io/travis/com/xumi1993/seispy-doc.post?label=Doc)](https://travis-ci.com/xumi1993/seispy-doc.post)
[![GitHub](https://img.shields.io/github/license/xumi1993/seispy)]()
[![](https://img.shields.io/github/last-commit/xumi1993/seispy)]()
[![](https://img.shields.io/github/commit-activity/m/xumi1993/seispy)]()
[![](https://img.shields.io/github/forks/xumi1993/seispy?style=social)]()

Python module of seismology and receiver functions

# Installation
## Dependencies
  * [Python]() >= 3.6
  * [ObsPy](http://docs.obspy.org) >= 1.1.0
  * [NumPy](http://www.numpy.org/) >= 1.16
  * [SciPy](http://www.scipy.org/) >= 1.2.0
  * [matplotlib](https://matplotlib.org/) >= 3.0.0
  * [PyQt5](https://www.riverbankcomputing.com/software/pyqt/)
  
## Installation
```Python
python setup.py install
```

# Inclusion
--------------
## Libraries
  * `seispy.distaz`: Calculate distance and azimuth (by [the lithospheric seismology program at USC](http://www.seis.sc.edu/software/distaz/)).<br />
  * `seispy.geo`: Tiny codes of geophysics.
  * `seispy.bootstrap`: Bootstrap confidence interval estimation (by [scikits-bootstrap](https://github.com/cgevans/scikits-bootstrap))
  * `seispy.decov`: Iterative time domain deconvolution method (Ligorria and Ammon's 1999 BSSA)

## Commands
 * `prf`: Calculate PRFs for a station.
 * `pickrf`: Pick PRFs after the calculation.
 * `plotrt`: Plot PRFs in R and T components order by back-azimuth.
 * `plotr`: Plot PRFs in R component order by back-azimuth.
 * `hk`: H-Kappa stacking.
 * `rf2depth`: Convert PRFs to depth axis.
 * `ccp_profile`: Stack PRFs along a profile with a CCP stacking method.

