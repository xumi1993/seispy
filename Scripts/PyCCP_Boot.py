#!/usr/bin/env python
import numpy as np
import scipy.io as sio
from scipy import interpolate
import os
import sys
import getopt
import seispy
import seispy.bootstrap as boot
import statsmodels.api as sm
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
# domperiod = float(config.get('para', 'domperiod'))
Profile_width = float(config.get('para', 'Profile_width'))
bin_radius = float(config.get('para', 'bin_radius'))
Stack_range = np.arange(300, 751)
smooth_para = 0.02
azi = seispy.distaz(lat1, lon1, lat2, lon2).baz
dis = seispy.distaz(lat1, lon1, lat2, lon2).delta
Profile_range = np.arange(0, seispy.geo.deg2km(dis), Profile_width)
Profile_lat = []
Profile_lon = []
for Profile_loca in Profile_range:
    (lat_loca, lon_loca) = seispy.geo.latlon_from(lat1, lon1, azi, seispy.geo.km2deg(Profile_loca))
#    Profile_lat = np.append(Profile_lat, [lat_loca], axis=1)
#    Profile_lon = np.append(Profile_lon, [lon_loca], axis=1)
    Profile_lat = np.append(Profile_lat, lat_loca)
    Profile_lon = np.append(Profile_lon, lon_loca)


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
    azi_sta = seispy.distaz(lat1, lon1, stalat_full[i], stalon_full[i]).baz
    dis_sta = seispy.distaz(lat1, lon1, stalat_full[i], stalon_full[i]).delta
    sta_dis = dis_sta * seispy.geo.sind(np.abs(azi - azi_sta))
    if sta_dis < seispy.geo.km2deg(bin_radius+310):
        print(RFdepth[0, i]['Station'][0], stalat_full[i], stalon_full[i])
        staidx = np.append(staidx, [i])

# stacking
Stack_data = np.zeros([Profile_range.shape[0],7], dtype=np.object)
#fid = open(stackfile, 'w+')
for i in range(Profile_range.shape[0]):
#    fid.write('>\n')
    print('calculate the RF stacks at the distance of '+str(Profile_range[i])+' km along the profile-------')
    mu_arr = np.zeros([Stack_range.shape[0],1])
    ci_arr = np.zeros([Stack_range.shape[0],2])
    Event_count_arr = np.zeros([Stack_range.shape[0],1])
    for j in range(Stack_range.shape[0]):
        print('No.'+str(i)+'/'+str(Profile_range.shape[0])+'---- at '+str(Stack_range[j])+'km depth.')
        Event_count = 0
        Amp_bin = []
#        if j == 0:
#            mu = np.nan
#            ci = np.array([np.nan, np.nan])
#            fid.write('%-8.3f %-8.3f %-6.3f %-6.3f %-8.5f %-8.5f %-8.5f %d\n' % (Profile_lat[i], Profile_lon[i], Profile_range[i], Stack_range[j], mu, ci[0], ci[1], Event_count))
#            mu_arr[j][0] = np.nan
#            ci_arr[j] = np.array([np.nan, np.nan])
#            Event_count_arr[j][0] = 0
#            continue
        for k in staidx:
            for l in range(RFdepth[0, k]['Piercelat'].shape[1]):
                # azi_event = distaz.distaz(Profile_lat[i], Profile_lon[i], RFdepth[0, k]['Piercelat'][l, 0], RFdepth[0, k]['Piercelon'][l, 0]).baz
                dis_event = seispy.distaz(Profile_lat[i], Profile_lon[i], RFdepth[0, k]['Piercelat'][2*Stack_range[j], l], RFdepth[0, k]['Piercelon'][2*Stack_range[j], l]).delta
                # pierce_interval = distaz.deg2km(dis_event) * np.abs(distaz.cosd(azi - azi_event))
                if dis_event < seispy.geo.km2deg(bin_radius):
                    Amp_bin = np.append(Amp_bin, RFdepth[0, k]['moveout_correct'][2*Stack_range[j], l])
                    Event_count = Event_count + 1
        if Event_count > 1:
            straptimes = int(Event_count*2)
            ci = boot.ci(Amp_bin, n_samples=straptimes, method='pi')
            mu = np.average(Amp_bin)
        else:
            mu = 0
            ci = np.zeros(2)
#        fid.write('%-8.3f %-8.3f %-6.3f %-6.3f %-8.5f %-8.5f %-8.5f %d\n' % (Profile_lat[i], Profile_lon[i], Profile_range[i], Stack_range[j], mu, ci[0], ci[1], Event_count))
        mu_arr[j][0] = mu
        ci_arr[j] = ci
        Event_count_arr[j][0] = Event_count
    mu_arr = sm.nonparametric.lowess(mu_arr[:, 0], Stack_range, frac=smooth_para)[:, 1]
    ci_arr[:, 0] = sm.nonparametric.lowess(ci_arr[:, 0], Stack_range, frac=smooth_para)[:, 1]
    ci_arr[:, 1] = sm.nonparametric.lowess(ci_arr[:, 1], Stack_range, frac=smooth_para)[:, 1]
    Stack_data[i][0] = Profile_lat[i]
    Stack_data[i][1] = Profile_lon[i]
    Stack_data[i][2] = Profile_range[i]
    Stack_data[i][3] = Stack_range
    Stack_data[i][4] = mu_arr
    Stack_data[i][5] = ci_arr
    Stack_data[i][6] = Event_count_arr
sio.savemat(stackfile,{'Stack_data':Stack_data})
