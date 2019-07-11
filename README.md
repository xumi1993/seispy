# seispy

[![pipeline status](https://git.nju.edu.cn/geophy/seispy/badges/master/pipeline.svg)](https://git.nju.edu.cn/geophy/seispy/commits/master)
[![LGPLv3](https://www.gnu.org/graphics/lgplv3-88x31.png)](https://www.gnu.org/licenses/lgpl.html)

Python module of seismology and  receiver functions

# Installation
## Dependencies
  * [Python]() >= 3.6
  * [ObsPy](http://docs.obspy.org) >= 1.1.0
  * [NumPy](http://www.numpy.org/) >= 1.14
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
  * `seispy.decov`: Iterative time domain deconv using Ligorria and Ammon's 1999 BSSA method

## Commands
 * `pickrf`: Pick PRFs after the calculation.
 * `plotrt`: Plot PRFs in R and T components order by back-azimuth.
 * `plotr`: Plot PRFs in R component order by back-azimuth.
 * `hk`: H-Kappa stacking.
 * `rf2depth`: Convert PRFs to depth axis.
 * `ccp_profile`: CCP stack PRFs along a profile.

