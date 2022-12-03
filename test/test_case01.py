import os
from seispy.rf import RF, read_catalog
from seispy.recalrf import ReRF
from seispy.hk import hksta
from seispy.hkpara import HKPara
from seispy.catalog import download_catalog
from subprocess import Popen
import pytest
from os.path import exists, join, basename
import glob
from obspy.io.sac import SACTrace
from obspy import UTCDateTime


def test_download():
    if exists('ex-prf.tar.gz'):
        pytest.skip('Data are downloaded.')
    s = 'wget https://osf.io/dxcfz/download -O ex-prf.tar.gz\n'
    s += 'tar -xzf ex-prf.tar.gz\n'
    proc = Popen(s, shell=True)
    proc.communicate()

def init_RF():
    rf = RF(cfg_file='ex-prf/rf.cfg')
    rf.para.phase = 'P'
    rf.para.datapath = 'ex-prf/Data.CB.NJ2'
    rf.para.rfpath = 'ex-prf/RFresult/CB.NJ2'
    rf.load_stainfo()
    rf.search_eq()
    rf.match_eq()
    rf.detrend()
    rf.filter()
    rf.cal_phase()
    rf.drop_eq_snr()
    rf.baz_correct()
    return rf

def gen_list(para):
    with open(os.path.join(para.rfpath, "CB.NJ2finallist.dat"), 'w+') as fid:
        files = sorted(glob.glob(join(para.rfpath, '*R.sac')))
        for fname in files:
            sac = SACTrace.read(fname)
            evname = basename(fname).split('_')[0]
            fid.write('%s %s %6.3f %6.3f %6.3f %6.3f %6.3f %8.7f %6.3f %6.3f\n' % (
                evname, 'P', sac.evla, sac.evlo, sac.evdp, sac.gcarc, sac.baz, sac.user0, sac.mag, sac.user1
            ))


def test_sub01():
    rf = init_RF()
    rf.para.decon_method = 'water'
    rf.para.criterion = 'crust'
    rf.para.rmsgate = None
    rf.rotate()
    rf.trim()
    rf.deconv()
    rf.saverf()


def test_sub02():
    rf = init_RF()
    rf.para.comp = 'lqt'
    rf.para.decon_method = 'water'
    rf.para.criterion = None
    rf.para.rmsgate = None
    rf.rotate()
    rf.trim()
    rf.deconv()
    rf.saverf()


def test_sub03():
    rf = init_RF()
    rf.para.decon_method = 'iter'
    rf.para.criterion = 'crust'
    rf.para.rmsgate = 0.2
    rf.rotate()
    rf.trim()
    rf.deconv()
    rf.saverf()
    gen_list(rf.para)


def test_sub04():
    rf = ReRF('ex-prf/RFresult/CB.NJ2/CB.NJ2finallist.dat',
              cfg_file='ex-prf/rf.cfg')
    rf.para.phase = 'P'
    rf.para.datapath = 'ex-prf/Data.CB.NJ2'
    rf.para.rfpath = 'ex-prf/RFresult/CB.NJ2_re'
    rf.para.decon_method = 'iter'
    rf.para.gauss = 2.5
    rf.para.criterion = None
    rf.para.rmsgate = None
    rf.load_stainfo()
    rf.match_eq()
    rf.detrend()
    rf.filter()
    rf.cal_phase()
    rf.rotate()
    rf.trim()
    rf.deconv()
    rf.saverf()
    rf.write_list()


def test_sub05():
    hkpara = HKPara()
    hkpara.rfpath = 'ex-prf/RFresult/CB.NJ2'
    hkpara.hkpath = './'
    hksta(hkpara, isplot=True, isdisplay=True)


def test_sub06():
    fname = 'evts.lst'
    b_time = UTCDateTime('20200101')
    e_time = UTCDateTime('20200201')
    download_catalog(fname, format='QUAKEML', starttime=b_time,
                     endtime=e_time, minmagnitude=5.5)
    eq_lst = read_catalog(fname, b_time, e_time, 0.,0.,)
    print(eq_lst)


if __name__ == '__main__':
    test_sub06()