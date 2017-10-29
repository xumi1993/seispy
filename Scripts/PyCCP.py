#!/usr/bin/env python

import numpy as np
import scipy.io as sio
from scipy import interpolate
import os
import sys
import getopt
import seispy
import seispy.bootstrap as boot
try:
    import configparser
    config = configparser.ConfigParser()
except:
    import ConfigParser
    config = ConfigParser.ConfigParser()

# get options of line location
try:
    opts, args = getopt.getopt(sys.argv[1:], "l:B")
except:
    print('Arguments are not found!')
    sys.exit(1)

isboot = False
for op, value in opts:
    if op == "-l":
        line = value
    elif op == "-B":
        isboot = True
    else:
        sys.exit(1)

# get parameters from cfg file
for o in sys.argv[1:]:
    if os.path.isfile(o):
        head = o
        config.read(head)
        break

# read parameters
lat1 = float(line.split('/')[0])
lon1 = float(line.split('/')[1])
lat2 = float(line.split('/')[2])
lon2 = float(line.split('/')[3])
depthdat = config.get('FileIO', 'depthdat')
stackfile = config.get('FileIO', 'stackfile')
Velmod = config.get('FileIO', 'Velmod')
stalist = config.get('FileIO', 'stalist')
domperiod = float(config.get('para', 'domperiod'))
Profile_width = float(config.get('para', 'Profile_width'))
# proj_width = float(config.get('para', 'proj_width'))
Stack_range = np.arange(1, 101)
azi = seispy.distaz(lat1, lon1, lat2, lon2).baz
dis = seispy.distaz(lat1, lon1, lat2, lon2).delta
Profile_range = np.arange(0, seispy.geo.deg2km(dis), Profile_width)
Profile_lat = []
Profile_lon = []
for Profile_loca in Profile_range:
    (lat_loca, lon_loca) = seispy.geo.latlon_from(lat1, lon1, azi, seispy.geo.km2deg(Profile_loca))
    #print(lat_loca, lon_loca)
    Profile_lat = np.append(Profile_lat, [lat_loca])
    Profile_lon = np.append(Profile_lon, [lon_loca])

# read depth data
depthmat = sio.loadmat(depthdat)
RFdepth = depthmat['YN_RFdepth']

# find stations beside the profile
stalat_all = RFdepth[0, 0::]['stalat']
stalon_all = RFdepth[0, 0::]['stalon']
stalat_full = []
stalon_full = []
staidx = []
for i in range(stalon_all.shape[0]):
    stalat_full = np.append(stalat_full, stalat_all[i][0])
    stalon_full = np.append(stalon_full, stalon_all[i][0])
STA = open(stalist, 'w+')
for i in range(stalon_all.shape[0]):
    (projstla, projstlo) = seispy.geo.geoproject(stalat_all[i][0], stalon_all[i][0],
                                                 lat1, lon1, lat2, lon2)
    sta_dis = seispy.distaz(projstla, projstlo, stalat_all[i][0], stalon_all[i][0]).delta
    if sta_dis < 0.5:
        print(RFdepth[0, i]['Station'][0], projstla, projstlo)
        STA.write('%s %-6.3f %-6.3f\n' % (RFdepth[0, i]['Station'][0], projstla, projstlo))
        staidx.append(i)
STA.close()

# Model information
YAxisRange = RFdepth[0, 1]['Depthrange'][0]
VelocityModel = np.loadtxt(Velmod)
Depths = VelocityModel[:, 0]
Vs = VelocityModel[:, 2]
Vs = interpolate.interp1d(Depths, Vs, kind='linear')(YAxisRange)

# stacking
fid = open(stackfile, 'w+')
for i in range(Profile_range.shape[0]):
    fid.write('>\n')
    print('calculate the RF stacks at the distance of '+str(Profile_range[i])+' km along the profile-------')
    dis_center = seispy.distaz(lat1, lon1, Profile_lat[i], Profile_lon[i]).degreesToKilometers()
    for j in range(Stack_range.shape[0]):
        bin_radius = np.sqrt(0.5*domperiod*Vs[2*Stack_range[j]]*Stack_range[j])
        Stack_RF = []
        Event_count = 0
        mu = 0
        for k in staidx:
            for l in np.arange(RFdepth[0, k]['Piercelat'].shape[1]):
                pier_lat = RFdepth[0, k]['Piercelat'][2*Stack_range[j], l]
                pier_lon = RFdepth[0, k]['Piercelon'][2*Stack_range[j], l]
                pier_azi = seispy.distaz(lat1, lon1, pier_lat, pier_lon).baz
                pier_dis = seispy.distaz(lat1, lon1, pier_lat, pier_lon).degreesToKilometers()
                dis_along = pier_dis * seispy.geo.cosd(azi - pier_azi)
                dis_proj = pier_dis * np.abs(seispy.geo.sind(azi - pier_azi))
                # (pro_lat, pro_lon) = seispy.geo.latlon_from(lat1, lon1, azi, seispy.geo.km2deg(dis_along))
                # dis_event = seispy.distaz(Profile_lat[i], Profile_lon[i], RFdepth[0, k]['Piercelat'][2*Stack_range[j], l], RFdepth[0, k]['Piercelon'][2*Stack_range[j], l]).degreesToKilometers()
                # pierce_interval = distaz.deg2km(dis_event) * np.abs(distaz.cosd(azi - azi_event))
                if np.abs(dis_along - Profile_range[i]) < bin_radius and dis_proj < seispy.geo.deg2km(bin_radius):
#                    print(pierce_interval, azi - azi_event, bin_radius)
#                    Stack_RF = Stack_RF + RFdepth[0, k]['moveout_correct'][2*Stack_range[j], l]
                    Stack_RF = np.append(Stack_RF, RFdepth[0, k]['moveout_correct'][2*Stack_range[j], l])
                    Event_count = Event_count + 1

        if isboot:
            if Event_count > 1:
                straptimes = int(Event_count * 2)
                ci = boot.ci(Stack_RF, n_samples=straptimes, method='pi')
                mu = np.average(Stack_RF)
            else:
                ci = np.zeros(2)
            if j == 0:
                mu = np.nan
                ci = np.array([np.nan, np.nan])
            fid.write('%-8.3f %-8.3f %-6.3f %-6.3f %-8.3f %-8.3f %-8.3f %d\n' %
                     (Profile_lat[i], Profile_lon[i], Profile_range[i], Stack_range[j], mu, ci[0], ci[1], Event_count))
        else:
            if Event_count > 0:
                mu = np.mean(Stack_RF)
            fid.write('%-8.3f %-8.3f %-6.3f %-6.3f %-8.3f %d\n' %
                    (Profile_lat[i], Profile_lon[i], Profile_range[i], Stack_range[j], mu, Event_count))
