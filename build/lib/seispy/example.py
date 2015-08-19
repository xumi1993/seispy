# -*- coding: utf-8 -*-
"""
Created on Wed Sep 10 22:43:21 2014

@author: xumj
"""

import decov
import numpy as np
import matplotlib.pyplot as plot

r=np.loadtxt('2012.045.06.17.01_RFdata_R.dat')
z=np.loadtxt('2012.045.06.17.01_RFdata_Z.dat')

dt=0.01
shift=10
itmax=400
mind=0.001
RFlength=13001
f0=2.0

RF=decov.decovit(r,z,dt,RFlength,shift,f0,itmax,mind)
Time=np.linspace(-shift, (RFlength-1)*dt-shift,RFlength)

plot.plot(Time,RF)
plot.plot([-2,30],[0,0])
plot.xlim(-2,30)
plot.show()
