#!/usr/bin/env python
import numpy as np
import scipy.io as sio
from scipy import interpolate
import os
import sys
import seispy
try:
    import configparser
    config = configparser.ConfigParser()
except:
    import ConfigParser
    config = ConfigParser.ConfigParser()


# get parameters from cfg file
for o in sys.argv[1:]:
    if os.path.isfile(o):
        head = o
        config.read(head)
        break

# read parameters
center_bin = config.get('FileIO', 'center_bin_local')
depthdat = config.get('FileIO', 'depthdat')
Velmod = config.get('FileIO', 'Velmod')
# domperiod = float(config.get('para', 'domperiod'))
bin_radius = float(config.get('para', 'bin_radius'))
Out_path = config.get('FileIO', 'outmat')
Stack_range = np.arange(0, 101)
center = np.loadtxt(center_bin)

# read depth data
depthmat = sio.loadmat(depthdat)
RFdepth = depthmat['YN_RFdepth']

# Model information
YAxisRange = RFdepth[0, 1]['Depthrange'][0]
VelocityModel = np.loadtxt(Velmod)
Depths = VelocityModel[:, 0]
Vs = VelocityModel[:, 2]
Vs = interpolate.interp1d(Depths, Vs, kind='linear')(YAxisRange)

# Stacking
Stack_data = np.zeros([center.shape[0],5], dtype=np.object)
for i in range(center.shape[0]):   
   print('Calculating %dth bin at %f/%f' % (i, center[i][0], center[i][1]))
   Stack_RF = np.zeros([Stack_range.shape[0], 1])
   Event_count = np.zeros([Stack_range.shape[0], 1])
   for j in range(Stack_range.shape[0]):
      # bin_radius = np.sqrt(0.5*domperiod*Vs[2*Stack_range[j] ]*Stack_range[j])
      for k in range(RFdepth.shape[1]):
         if seispy.distaz(center[i][0], center[i][1], RFdepth[0, k]['stalat'][0][0], RFdepth[0, k]['stalon'][0][0]).delta>=2:
            continue
         for l in range(RFdepth[0, k]['Piercelat'].shape[1]):
            if seispy.distaz(center[i][0], center[i][1], RFdepth[0, k]['Piercelat'][2*Stack_range[j], l], RFdepth[0, k]['Piercelon'][2*Stack_range[j], l]).degreesToKilometers()<bin_radius:
               Stack_RF[j][0] += RFdepth[0, k]['moveout_correct'][2*Stack_range[j], l]
               Event_count[j][0] += 1
      if Event_count[j][0] > 0:
         Stack_RF[j][0] /= Event_count[j][0]
   Stack_data[i][0] = center[i][0]
   Stack_data[i][1] = center[i][1]
   Stack_data[i][2] = Stack_range
   Stack_data[i][3] = Stack_RF
   Stack_data[i][4] = Event_count
sio.savemat(Out_path,{'Stack_data':Stack_data})

