#!/usr/bin/env python

def Usage():
   print('./initgrid.py -Rlon1/lon2/lat1/lat2 -dlonval/latval')

import os
import numpy as np
import getopt
import sys

try:
   opts, args = getopt.getopt(sys.argv[1:], "R:d:")
except:
    print('Arguments are not found!')
    Usage()
    sys.exit(1)

if opts == []:
   Usage()
   sys.exit(1)

for op, value in opts:
   if op == "-R":
      RridRange = value
   elif op == "-d":
      Interval = value
   else:
      Usage()
      sys.exit(1)

lon1 = float(RridRange.split('/')[0])
lon2 = float(RridRange.split('/')[1])
lat1 = float(RridRange.split('/')[2])
lat2 = float(RridRange.split('/')[3])

LonVal = float(Interval.split('/')[0])
LatVal = float(Interval.split('/')[1])

LonRange = np.arange(lon1,lon2+LonVal,LonVal)
LatRange = np.arange(lat1,lat2+LatVal,LatVal)

for lon in LonRange:
   for lat in LatRange:
      print(lat, lon)

