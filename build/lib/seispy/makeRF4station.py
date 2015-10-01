#!/usr/bin/env python

def Usage():
   print('./makeRF4station -Sstation/slat/slon -Yyear1/month1/day1/year2/month2/day2 -Ccomp -T+N+E para.cfg')

import numpy as np
import obspy
import re
import os
import glob
import getopt
import sys
import datetime
import seispy
from obspy.taup import TauPyModel
from scipy.signal import detrend
from obspy.signal.filter import bandpass
from scipy.signal import resample
try:
    import configparser
    config = configparser.ConfigParser()
except:
    import ConfigParser
    config = ConfigParser.ConfigParser()

#################################
# ------- Get arguments --------#
#################################
try:
   opts, args = getopt.getopt(sys.argv[1:], "S:Y:C:T:")
except:
    print('Arguments are not found!')
    Usage()
    sys.exit(1)
if opts == []:
    Usage()
    sys.exit(1)

change_E = 0
change_N = 0
for op, value in opts:
    if op == "-S":
       station = value
    elif op == "-Y":
       date_str = value
    elif op == "-C":
        comp = value
    elif op == "-T":
        change_chan = value
        change_E = change_chan.find("E")
        change_N = change_chan.find("N")
    else:
       Usage()
       sys.exit(1)

for o in sys.argv[1:]:
    if os.path.isfile(o):
        head = o
        config.read(head)
        break

#################################
# --------- Get Para -----------#
#################################
model = TauPyModel(model="iasp91")
staname = station.split('/')[0]
stalat = float(station.split('/')[1])
stalon = float(station.split('/')[2])
yearrange1 = int(date_str.split('/')[0])
monthrange1 = int(date_str.split('/')[1])
dayrange1 = int(date_str.split('/')[2])
yearrange2 = int(date_str.split('/')[3])
monthrange2 = int(date_str.split('/')[4])
dayrange2 = int(date_str.split('/')[5])
daterange1 = datetime.datetime(yearrange1, monthrange1, dayrange1)
daterange2 = datetime.datetime(yearrange2, monthrange2, dayrange2)
data_path = config.get('path', 'data_path')
out_path = config.get('path', 'out_path')
out_path = os.path.join(out_path, staname)
RF_path = config.get('path', 'RF_path')
RF_path = os.path.join(RF_path, staname)
evt_list = config.get('path', 'evt_list')
image_path = config.get('path', 'image_path')
if not os.path.exists(image_path):
    os.makedirs(image_path)
gate_mw = float(config.get('para', 'gate_mw'))
gate_dis1 = float(config.get('para', 'gate_dis1'))
gate_dis2 = float(config.get('para', 'gate_dis2'))
time_before = float(config.get('para', 'time_before'))
time_after = float(config.get('para', 'time_after'))
searchchdt = float(config.get('para', 'searchchdt'))
offset = float(config.get('para', 'offset'))
gate_noise = float(config.get('para', 'gate_noise'))
gauss = float(config.get('para', 'gauss'))
freqmin = float(config.get('para', 'freqmin'))
freqmax = float(config.get('para', 'freqmax'))
sampling = float(config.get('para', 'sampling'))
newlength = ((time_before+time_after)/sampling)+1
image_path = os.path.join(image_path, staname+'_R_gauss'+str(gauss)+'.ps')
fid_evtlist = open(evt_list, 'r')

#################################
# --------- Search eq ----------#
#################################
eq_lst = []
for evt in fid_evtlist.readlines():
    evt = evt.strip('\n')
    year_evt = int(evt.split()[0])
    mon_evt = int(evt.split()[1])
    day_evt = int(evt.split()[2])
    hour_evt = int(evt.split()[4])
    min_evt = int(evt.split()[5])
    sec_evt = int(evt.split()[6])
    lat_evt = float(evt.split()[7])
    lon_evt = float(evt.split()[8])
    dis_evt = seispy.distaz(stalat, stalon, lat_evt, lon_evt).delta
    bazi_evt = seispy.distaz(stalat, stalon, lat_evt, lon_evt).getBaz()
    dep_evt = float(evt.split()[9])
    mw_evt = float(evt.split()[10])
    date_evt = datetime.datetime(year_evt, mon_evt, day_evt, hour_evt, min_evt, sec_evt)
    if daterange1 <= date_evt <= daterange2 and gate_dis1 <= dis_evt <= gate_dis2 and mw_evt > gate_mw:
        matchedeq = [date_evt, lat_evt, lon_evt, dep_evt, dis_evt, bazi_evt, mw_evt]
        eq_lst.append(matchedeq)

#################################
# -------- Assign Files --------#
#################################
eq = []
for nowevt in eq_lst:
    find_result = []
    for sac in glob.glob(os.path.join(data_path, staname, '*.'+comp+'.*.[Ss][Aa][Cc]')):
        sacname = os.path.basename(sac)
        date_name = re.search("\d{4}\D\d{3}\D\d{2}\D\d{2}\D\d{2}", sacname).group()
        date_name_sp = re.split('\D', date_name)
        year_sac = int(date_name_sp[0])
        jday_sac = int(date_name_sp[1])
        foo = datetime.datetime(year_sac, 1, 1) + datetime.timedelta(jday_sac - 1)
        mon_sac = foo.month
        day_sac = foo.day
        hour_sac = int(date_name_sp[2])
        min_sac = int(date_name_sp[3])
        sec_sac = int(date_name_sp[4])
        if sec_sac == 60:
            sec_sac = 0
            date_sac = datetime.datetime(year_sac, mon_sac, day_sac, hour_sac, min_sac, sec_sac) + datetime.timedelta(seconds=60)
        else:
            date_sac = datetime.datetime(year_sac, mon_sac, day_sac, hour_sac, min_sac, sec_sac)
        # if nowevt[0] + datetime.timedelta(seconds=offset) - datetime.timedelta(seconds=searchchdt) <= date_sac <= nowevt[0] + datetime.timedelta(seconds=offset) + datetime.timedelta(seconds=searchchdt):
        if date_sac - datetime.timedelta(seconds=offset) - datetime.timedelta(seconds=searchchdt) <= nowevt[0] <= date_sac - datetime.timedelta(seconds=offset) + datetime.timedelta(seconds=searchchdt):
            print(date_sac, nowevt[0])
            ocr = (date_sac - nowevt[0])
            ocr = ocr.total_seconds()
            find_result.append([nowevt[0], date_name, nowevt[1], nowevt[2], nowevt[3], nowevt[4], nowevt[5], nowevt[6], ocr])
    if len(find_result) == 1:
        eq.append(find_result[0])
print(len(eq))

#################################
# -------- Process RFs ---------#
#################################
if not os.path.exists(RF_path):
    os.makedirs(RF_path)
fid_filallist = open(os.path.join(RF_path, staname+'finallist.dat'), 'w+')
for thiseq in eq:
    date_name = thiseq[1]
    bazi = thiseq[6]
    dis = thiseq[5]
    dep = thiseq[4]
    O = thiseq[-1]
    try:
        this_seis = obspy.core.read(os.path.join(data_path, staname, '*'+date_name+'*.[Ss][Aa][Cc]'))
    except:
        continue
# --------- Calculate dt -----------#
    try:
        dt_E = np.mean(np.diff(this_seis[0].times()))
        dt_N = np.mean(np.diff(this_seis[1].times()))
        dt_Z = np.mean(np.diff(this_seis[2].times()))
    except:
        continue
    if not (dt_E == dt_N == dt_Z and len(this_seis[0].data) == len(this_seis[1].data) == len(this_seis[2].data)):
        continue
    dt = dt_E
# -------- Detrend -----------------#
    if change_E:
        this_seis[0].data = this_seis[0].data * -1
    if change_N:
        this_seis[1].data = this_seis[1].data * -1
    this_seis[0].data = detrend(this_seis[0].data, type='linear')
    this_seis[0].data = detrend(this_seis[0].data, type='constant')
    this_seis[1].data = detrend(this_seis[1].data, type='linear')
    this_seis[1].data = detrend(this_seis[1].data, type='constant')
    this_seis[2].data = detrend(this_seis[2].data, type='linear')
    this_seis[2].data = detrend(this_seis[2].data, type='constant')
# ------- Rotate ENZ to RTZ ---------#
    if comp == 1 or comp == 2 or comp == 3:
        E = this_seis[2].data
        N = this_seis[1].data
        Z = this_seis[0].data
    else:
        E = this_seis[0].data
        N = this_seis[1].data
        Z = this_seis[2].data
    (this_seis[0].data, this_seis[1].data, this_seis[2].data) = seispy.geo.rotateSeisENZtoTRZ(E, N, Z, bazi)
    this_seis[0].stats.channel = "R"
    this_seis[1].stats.channel = "T"
    this_seis[2].stats.channel = "Z"
# ------- Bandpass Filter ------------#
    R_filter = bandpass(this_seis[0].data, freqmin, freqmax, 1/dt, corners=3, zerophase=True)
    T_filter = bandpass(this_seis[1].data, freqmin, freqmax, 1/dt, corners=3, zerophase=True)
    Z_filter = bandpass(this_seis[2].data, freqmin, freqmax, 1/dt, corners=3, zerophase=True)
# ------- Calculate P arrival time ---#
    arrivals = model.get_travel_times(source_depth_in_km=dep, distance_in_degree=dis, phase_list=["P"])
    ttime_P = arrivals[0].time - O
    rayp = seispy.geo.srad2skm(arrivals[0].ray_param)
    # this_seis[0].stats.sac.t1 = ttime_P
    # this_seis[1].stats.sac.t1 = ttime_P
    # this_seis[2].stats.sac.t1 = ttime_P
# ----------- filter by snr ----------#
    snrbegin = np.floor((ttime_P-50)/dt)
    snrend = np.floor((ttime_P+50)/dt)
    snro = np.floor(ttime_P/dt)
    snr_R = seispy.geo.snr(R_filter[snro:snrend], R_filter[snrbegin:snro])
    snr_T = seispy.geo.snr(T_filter[snro:snrend], T_filter[snrbegin:snro])
    snr_Z = seispy.geo.snr(Z_filter[snro:snrend], Z_filter[snrbegin:snro])
    print(date_name, snr_R, snr_Z, ttime_P, O)
    if snr_R > gate_noise and snr_Z > gate_noise:
        extbegin = np.floor((ttime_P-time_before)/dt)
        extend = np.floor((ttime_P+time_after)/dt)
        R = R_filter[extbegin:extend+1]
        T = T_filter[extbegin:extend+1]
        Z = Z_filter[extbegin:extend+1]
        RFlength = R.shape[0]
# --- calculate receiver functions ----#
        (RF, RMS, it) = seispy.decov.decovit(R, Z, dt, RFlength, time_before, gauss, 400, 0.001)
        max_deep = np.max(np.abs(RF[np.floor((25+time_before)/dt):-1]))
        time_P1 = np.floor((-2+time_before)/dt)
        time_P2 = np.floor((2+time_before)/dt)
        max_P = np.max(RF[time_P1:time_P2])
        if RMS[-1] < 0.2 and max_deep < max_P*0.3 and max_P == np.max(np.abs(RF)) and max_P < 1:
            print("----"+staname+"-----"+date_name+"-----")
            if RFlength != newlength:
                RF = resample(RF, newlength)
            this_seis[0].data = this_seis[0].data[extbegin:extend+1]
            this_seis[1].data = this_seis[1].data[extbegin:extend+1]
            this_seis[2].data = this_seis[2].data[extbegin:extend+1]
            timeaxis = np.linspace(-time_before, time_after, RFlength)
            for i in range(3):
                this_seis[i].stats.sac.evla = thiseq[2]
                this_seis[i].stats.sac.evlo = thiseq[3]
                this_seis[i].stats.sac.evdp = dep
                this_seis[i].stats.sac.baz = bazi
                this_seis[i].stats.sac.gcarc = dis
                this_seis[i].stats.sac.user0 = rayp
                this_seis[i].stats.sac.t1 = time_before
            if not os.path.exists(out_path):
                os.makedirs(out_path)
            this_seis[0].write(os.path.join(out_path, date_name+'.'+staname+'.R.SAC'), 'SAC')
            this_seis[1].write(os.path.join(out_path, date_name+'.'+staname+'.T.SAC'), 'SAC')
            this_seis[2].write(os.path.join(out_path, date_name+'.'+staname+'.Z.SAC'), 'SAC')
            np.savetxt(os.path.join(RF_path, date_name+'_P_R.dat'), RF, fmt='%10.8f')
            fid_filallist.write('%s %s %6.3f %6.3f %6.3f %6.3f %6.3f %8.7f %6.3f %6.3f\n' % (date_name, 'P', thiseq[2], thiseq[3], dep, dis, bazi, rayp, thiseq[7], gauss))
fid_filallist.close()
# os.popen('cp '+os.path.join(RF_path, staname+'finallist.dat')+' '+out_path)
# os.popen('python /Users/xumj/Researches/PyCCP/PlotR.py -S%s -I%s -O%s -N%d -A80 -T%5.2f/%5.2f' % (staname, RF_path, image_path, newlength, time_before, time_after))
