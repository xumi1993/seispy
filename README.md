# seispy

[![pipeline status](https://git.nju.edu.cn/geophy/seispy/badges/master/pipeline.svg)](https://git.nju.edu.cn/geophy/seispy/commits/master)

Python module of seismology and  receiver functions

# Installation
## Dependencies
  * [obspy](http://docs.obspy.org)<br />
  * [GMT 5.x](http://gmt.soest.hawaii.edu)<br />
  
## Installation
```Python
python setup.py install
```
Add the ```seispy/Scripts``` to environment	variables.
# Inclusion
--------------
  * seispy.distaz: Calculate distance and azimuth (by [the lithospheric seismology program at USC](http://www.seis.sc.edu/software/distaz/)).<br />
  * seispy.geo: Tiny codes of geophysics.
  * seispy.bootstrap: Bootstrap confidence interval estimation (by [scikits-bootstrap](https://github.com/cgevans/scikits-bootstrap))
  * seispy.decov: Iterative time domain deconv using Ligorria and Ammon's 1999 BSSA method
