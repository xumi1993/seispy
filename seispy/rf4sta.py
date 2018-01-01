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


def searcheq(stalat, stalon, daterange1, daterange2, gate_dis1, gate_dis2, gate_mw, evt_list):
    eq_lst = []
    fid_evtlist = open(evt_list)
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
    return eq_lst


def assignFiles(files, st, eq_lst, tolerance):
    eq = []
    for nowevt in eq_lst:
        find_result = []
        for filename, tr in zip(files, st):
            sac_otime = tr.stats.starttime + datetime.timedelta(seconds=int(tr.stats.sac.o))
            if sac_otime - datetime.timedelta(seconds=tolerance) <= nowevt[0] <= sac_otime + datetime.timedelta(seconds=tolerance):
                print(sac_otime, nowevt[0])
                find_result.append(
                    [nowevt[0], filename, nowevt[1], nowevt[2], nowevt[3], nowevt[4], nowevt[5], nowevt[6]])
        if len(find_result) == 1:
            eq.append(find_result[0])
    return eq



class starf():
    def __init__(self, cfg_file, staname, suffix, daterange1, daterange2):
        self.model = TauPyModel(model="iasp91")
        self.daterange1 = daterange1
        self.daterange2 = daterange2
        config.read(cfg_file)
        self.data_path = os.path.join(config.get('path', 'data_path'), staname)
        self.out_path = os.path.join(config.get('path', 'out_path'), staname)
        self.RF_path = os.path.join(config.get('path', 'RF_path'), staname)
        self.evt_list = config.get('path', 'evt_list')
        self.image_path = config.get('path', 'image_path')
        self.gate_mw = config.getfloat('para', 'gate_mw')
        self.gate_dis1 = config.getfloat('para', 'gate_dis1')
        self.gate_dis2 = config.getfloat('para', 'gate_dis2')
        self.time_before = config.getfloat('para', 'time_before')
        self.time_after = config.getfloat('para', 'time_after')
        self.tolerance = config.getfloat('para', 'tolerance')
        self.offset = config.getfloat('para', 'offset')
        self.gate_noise = config.getfloat('para', 'gate_noise')
        self.gauss = config.getfloat('para', 'gauss')
        self.freqmin = config.getfloat('para', 'freqmin')
        self.freqmax = config.getfloat('para', 'freqmax')
        self.sampling = config.getfloat('para', 'sampling')
        self.newlength = ((self.time_before + self.time_after) / self.sampling) + 1

        self.zfiles = glob.glob(os.path.join(self.data_path, suffix))
        self.st = obspy.Stream()
        for zfile_name in self.zfiles:
            self.st.append(obspy.read(zfile_name)[0])
        print(len(self.st))
        self.stalat = self.st[0].stats.sac.stla
        self.stalon = self.st[0].stats.sac.stlo
        self.eq_lst = searcheq(self.stalat, self.stalon, self.daterange1, self.daterange2,
                               self.gate_dis1, self.gate_dis2, self.gate_mw, self.evt_list)
        self.eq = assignFiles(self.zfiles, self.st, self.eq_lst, self.tolerance)



if __name__ == '__main__':
    cfg_file = '/home/xumj/Codes/seispy/Scripts/paraRF.cfg'
    staname = 'KKD02'
    suffix = r'*.BHZ'
    datereange1 = obspy.UTCDateTime('20060101')
    datereange2 = obspy.UTCDateTime('20100101')
    rf = starf(cfg_file, staname, suffix, datereange1, datereange2)