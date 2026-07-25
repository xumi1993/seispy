import glob
import inspect
import os
from os.path import basename, join
from pathlib import Path
from shutil import copyfileobj
from tarfile import open as open_tar
from urllib.request import Request, urlopen

import pytest
from obspy import UTCDateTime
from obspy.io.sac import SACTrace

from seispy.catalog import download_catalog
from seispy.hk import hksta
from seispy.hkpara import HKPara
from seispy.recalrf import ReRF
from seispy.rf import RF, read_catalog


EX_PRF_URL = 'https://osf.io/dxcfz/download'
EX_PRF_ARCHIVE = Path('ex-prf.tar.gz')
EX_PRF_CONFIG = Path('ex-prf/rf.cfg')


def test_download():
    if EX_PRF_CONFIG.is_file():
        pytest.skip('Data are downloaded.')

    request = Request(EX_PRF_URL, headers={'User-Agent': 'seispy-ci'})
    with urlopen(request, timeout=120) as response:
        with EX_PRF_ARCHIVE.open('wb') as archive_file:
            copyfileobj(response, archive_file)

    with open_tar(EX_PRF_ARCHIVE, 'r:gz') as archive:
        extract_options = {}
        if 'filter' in inspect.signature(archive.extractall).parameters:
            extract_options['filter'] = 'data'
        archive.extractall('.', **extract_options)

    assert EX_PRF_CONFIG.is_file(), (
        'The ex-prf archive was downloaded but did not contain '
        f'{EX_PRF_CONFIG}.'
    )

def init_RF():
    rf = RF(cfg_file='ex-prf/rf.cfg')
    rf.para.phase = 'P'
    rf.para.datapath = 'ex-prf/Data.CB.NJ2'
    rf.para.rfpath = 'ex-prf/RFresult/CB.NJ2'
    rf.para.cata_server = 'USGS'
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
