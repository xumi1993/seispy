from seispy.rfcorrect import RFStation
import pytest
from os.path import exists
from subprocess import Popen
from seispy.plotR import plotr
from seispy.plotRT import plotrt


def test_download():
    if exists('ex-rfani.tar.gz'):
        pytest.skip('Data are downloaded.')
    s = 'wget https://osf.io/4hk6d/download -O ex-rfani.tar.gz\n'
    s += 'tar -xzf ex-rfani.tar.gz\n'
    proc = Popen(s, shell=True)
    proc.communicate()


def test_sub01():
    rfs = RFStation('ex-rfani/SC.LTA')
    rfs.jointani(3, 8, tlen=3.5, stack_baz_val=10, weight=[0.4, 0.4, 0.2])


def test_sub02():
    rfs = RFStation('ex-rfani/SC.LTA')
    rfs.moveoutcorrect()
    rfs.psrf2depth()
    rfs.psrf_1D_raytracing()
    rfs.slantstack()


def test_sub03():
    rfs = RFStation('ex-rfani/SC.LTA')
    rfs.resample(0.1)
    rfs.harmonic(-2, 12)
    rfs.harmo.write_constant()
    rfs.harmo.plot()


def test_sub04():
    rfs = RFStation('ex-rfani/SC.LTA')
    rfs.resample(0.1)
    plotr(rfs, outpath='./', xlim=[-2, 80], key='bazi', enf=6, format='pdf')
    plotrt(rfs, enf=3, out_path='./', key='bazi', outformat='g', xmax=30)