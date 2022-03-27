from seispy.rfcorrect import RFStation
import pytest
from os.path import exists
from subprocess import Popen


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
