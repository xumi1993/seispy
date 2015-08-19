#!/usr/bin/env python
import numpy as np
import scipy.io as sio
from scipy import interpolate
import os
import sys
import getopt
import distaz
import scikits.bootstrap as boot
try:
    import configparser
    config = configparser.ConfigParser()
except:
    import ConfigParser
    config = ConfigParser.ConfigParser()

# get options of line location
try:
    opts, args = getopt.getopt(sys.argv[1:], "l:")
except:
    print('Arguments are not found!')
    sys.exit(1)

for op, value in opts:
    if op == "-l":
        line = value
    else:
        sys.exit(1)

# get parameters from cfg file
for o in sys.argv[1:]:
    if os.path.isfile(o):
        head = o
        config.read(head)
        break

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
bin_radius = float(config.get('para', 'bin_radius'))
Stack_range = np.arange(350, 751)
azi = distaz.distaz(lat1, lon1, lat2, lon2).baz
dis = distaz.distaz(lat1, lon1, lat2, lon2).delta
Profile_range = np.arange(0, distaz.deg2km(dis), Profile_width)
Profile_lat = []
Profile_lon = []
for Profile_loca in Profile_range:
    (lat_loca, lon_loca) = distaz.latlon_from(lat1, lon1, azi, distaz.km2deg(Profile_loca))
    Profile_lat = np.append(Profile_lat, [lat_loca], axis=1)
    Profile_lon = np.append(Profile_lon, [lon_loca], axis=1)

# ----- Read depth .mat file -----#
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
for i in range(stalon_all.shape[0]):
    azi_sta = distaz.distaz(lat1, lon1, stalat_full[i], stalon_full[i]).baz
    dis_sta = distaz.distaz(lat1, lon1, stalat_full[i], stalon_full[i]).delta
    sta_dis = dis_sta * distaz.sind(np.abs(azi - azi_sta))
    if sta_dis < distaz.km2deg(bin_radius+310):
        print(RFdepth[0, i]['Station'][0], stalat_full[i], stalon_full[i])
        staidx = np.append(staidx, [i])

# stacking
fid = open(stackfile, 'w+')
for i in range(Profile_range.shape[0]):
    fid.write('>\n')
    print('calculate the RF stacks at the distance of '+str(Profile_range[i])+' km along the profile-------')
    for j in range(Stack_range.shape[0]):
        print('No.'+str(i)+'/'+str(Profile_range.shape[0])+'---- at '+str(Stack_range[j])+'km depth.')
        Event_count = 0
        Amp_bin = []
        if j == 0:
            mu = np.nan
            ci = np.array([np.nan, np.nan])
            fid.write('%-8.3f %-8.3f %-6.3f %-6.3f %-8.5f %-8.5f %-8.5f %d\n' % (Profile_lat[i], Profile_lon[i], Profile_range[i], Stack_range[j], mu, ci[0], ci[1], Event_count))
            continue
        for k in staidx:
            for l in range(RFdepth[0, k]['Piercelat'].shape[1]):
                # azi_event = distaz.distaz(Profile_lat[i], Profile_lon[i], RFdepth[0, k]['Piercelat'][l, 0], RFdepth[0, k]['Piercelon'][l, 0]).baz
                dis_event = distaz.distaz(Profile_lat[i], Profile_lon[i], RFdepth[0, k]['Piercelat'][2*Stack_range[j] + 1, l], RFdepth[0, k]['Piercelon'][2*Stack_range[j] + 1, l]).delta
                # pierce_interval = distaz.deg2km(dis_event) * np.abs(distaz.cosd(azi - azi_event))
                if dis_event < distaz.km2deg(bin_radius):
                    Amp_bin = np.append(Amp_bin, RFdepth[0, k]['moveout_correct'][2*Stack_range[j] + 1, l])
                    Event_count = Event_count + 1
        if Event_count > 1:
            straptimes = np.round(Event_count*1.5)
            print(Event_count,straptimes)
            ci = boot.ci(Amp_bin, n_samples=straptimes, method='pi')
            mu = np.average(Amp_bin)
        else:
            mu = 0
            ci = np.zeros(2)
        fid.write('%-8.3f %-8.3f %-6.3f %-6.3f %-8.5f %-8.5f %-8.5f %d\n' % (Profile_lat[i], Profile_lon[i], Profile_range[i], Stack_range[j], mu, ci[0], ci[1], Event_count))