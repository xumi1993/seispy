import os
from seispy.rf import RF
from obspy.signal.rotate import rotate2zne
from obspy import UTCDateTime
from obspy.io.sac import SACTrace
import numpy as np
#import pytest
import glob
import matplotlib.pyplot as plt

def main():
    #Set the parameters for RF calculation
    rf = RF()
    rf.para.data_server = 'IRIS'
    rf.para.cata_server = 'IRIS'
    # collected data at the following stations, do not forget to rotate to ZNE first and use correct location code::
    # G.SANVU: 20111101-20250101, location='00',
    # G.PVC: 19940601-20040401, location='--', not NE component
    # G.DZM:20030901-20250101, location='00'
    # G.NOUC: 19880321-20250101 (near G.DZM), location='00' or '--', so do not specify location
    # IU.FUNA: 20040121-20250101, location='00', not NE component
    # G.FUTU: 20160626-20250101, location='00'
    # II.MSVF: 19940524-20250101, location='00', not NE component
    # IU.AFI: 20000101-20250101, location='00' , not NE component (but close)
    # AU.NIUE: 20070621-20250101,location='00' or '--', so do not specify location
    # IU.RAO: 20040718-20250101, location='00'
    rf.para.stainfo.network = 'IU'
    rf.para.stainfo.station = 'AFI'
    rf.para.stainfo.location = '00'
    #rf.para.stainfo.location = '--'
    rf.para.stainfo.channel = 'BH?'
    rf.para.date_begin = UTCDateTime('20000101')
    rf.para.date_end = UTCDateTime('20250101')

    rf.para.datapath = './Data/{}.{}'.format(rf.para.stainfo.network, rf.para.stainfo.station)
    rf.para.use_remote_data = True
    rf.para.ref_comp ='BHZ'
    rf.para.phase = 'P'
    rf.para.noisegate = 5   # SNR criterion, SNR=10log10(Amp_signal/Amp_noise),
    # where Amp_signal is the RMS amplitude of the signal and Amp_noise is the RMS amplitude of the noise. 5 corresponds to amplitude ratio of 3.16.
    rf.para.magmin = 5.5
    rf.para.dismin = 30
    rf.para.dismax = 90
    rf.para.time_before = 10
    rf.para.time_after = 120
    rf.para.gauss = [0.5, 1.0, 1.5] ##RF with different Gauss factor will be calculated simultaneously.
    rf.para.rmsgate = 0.4
    rf.para.freqmin = 0.02
    rf.para.freqmax = 2.0
    rf.para.comp = 'RTZ'
    rf.para.n_proc = 8
    print(rf.para)

    #load station information and search events
    rf.load_stainfo()
    rf.search_eq(catalog='NEIC PDE')
    #print(rf.eq_lst) ##The matched event lists are listed.
    #Match catalog and fetch seismic data
    rf.match_eq()

    # preprocessing, quality control, and save raw data
    rf.detrend()  # this step contains rotating non-traditional component to ZNE component
    rf.filter()
    rf.cal_phase() # Get arrival time, ray parameter and incident angle from TauP model
    rf.rotate()
    # quality control, save, trim
    rf.drop_eq_snr()   # apply SNR criterion to drop low SNR events
    rf.save_raw_data()
    rf.trim()

    # calculate and save RFs
    rf.deconv()
    for ff in rf.para.gauss:
        rf.para.rfpath = './RFresult/F{:.1f}/{}.{}'.format(ff, rf.para.stainfo.network, rf.para.stainfo.station)
        rf.saverf(ff)

if __name__ == '__main__':
    main()


