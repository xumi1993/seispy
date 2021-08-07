# <img src="https://user-images.githubusercontent.com/7437523/128596331-dc5c5e40-93e1-4d9e-b92d-9c53fe51145a.png" width="500"/> 

[![pipeline status](https://img.shields.io/travis/com/xumi1993/seispy)](https://travis-ci.com/xumi1993/seispy)
[![Build Status](https://img.shields.io/travis/com/xumi1993/seispy-doc.post?label=doc)](https://seispy.xumijian.me)
[![PyPI](https://img.shields.io/pypi/v/python-seispy)](https://pypi.org/project/python-seispy/)
[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/python-seispy)](https://pypi.org/project/python-seispy/)
[![GitHub](https://img.shields.io/github/license/xumi1993/seispy)]()
[![GitHub tag (latest by date)](https://img.shields.io/github/v/tag/xumi1993/seispy)]()
[![](https://img.shields.io/github/last-commit/xumi1993/seispy)]()
[![](https://img.shields.io/github/commit-activity/m/xumi1993/seispy)]()
[![GitHub code size in bytes](https://img.shields.io/github/languages/code-size/xumi1993/seispy)]()
[![GitHub repo size](https://img.shields.io/github/repo-size/xumi1993/seispy)]()
[![DOI](https://zenodo.org/badge/41006349.svg)](https://zenodo.org/badge/latestdoi/41006349)

[![GitHub stars](https://img.shields.io/github/stars/xumi1993/seispy?style=social)]()
[![](https://img.shields.io/github/forks/xumi1993/seispy?style=social)]()


Seispy is a Python module for processing seismological data and calculating Receiver Functions. The advanced functions are available to improve the Obspy.


# Installation
## Dependencies
  * [Python]() >= 3.6
  * [ObsPy](http://docs.obspy.org) >= 1.1.0
  * [NumPy](http://www.numpy.org/) >= 1.16
  * [SciPy](http://www.scipy.org/) >= 1.2.0
  * [matplotlib](https://matplotlib.org/) >= 3.0.0
  * [PyQt5](https://www.riverbankcomputing.com/software/pyqt/) >= 5.12.0
  
## Installation
```
git clone https://github.com/xumi1993/seispy.git
python setup.py install
```

# Inclusion
## Libraries
  * `seispy.distaz`: Calculate distance and azimuth (by [the lithospheric seismology program at USC](http://www.seis.sc.edu/software/distaz/)).<br />
  * `seispy.geo`: Tiny codes of geophysics.
  * `seispy.bootstrap`: Bootstrap confidence interval estimation (by [scikits-bootstrap](https://github.com/cgevans/scikits-bootstrap))
  * `seispy.decon`: Iterative time domain deconvolution method (Ligorria and Ammon's 1999 BSSA)
  * `seispy.rfcorrect`: Subsequent process of PRFs includeing moveout correct and time to depth conversion (1D and 3D) (see [Mijian Xu et al., 2018 EPSL](https://www.sciencedirect.com/science/article/pii/S0012821X17306921?via%3Dihub))
  * `seispy.ccp`: CCP stacking along a profile.


## Commands
### Receiver Functions
 * `prf`: Calculate PRFs for a station.
 * `pickrf`: Pick PRFs with virtual quality control after the calculation.
 * `plotrt`: Plot PRFs in R and T components order by back-azimuth.
 * `plotr`: Plot PRFs in R component order by back-azimuth.
 * `hk`: H-Kappa stacking.
 * `rf2depth`: Convert PRFs to depth axis.
 * `ccp_profile`: Stack PRFs along a profile with a CCP stacking method.

### Others
 * `ndk2dat`: Convert the GCMT catalog file ("ndk" format) to the list file for the `prf` command.
 * `updatecatalog`: Automatically update the GCMT catalog.
 * `setpar`: Set up the values in configure files.

