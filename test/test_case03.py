from seispy.ccppara import ccppara
from seispy.ccpprofile import CCPProfile
from seispy.rf2depth_makedata import makedata
from os.path import exists
from subprocess import Popen
import pytest
import os

def test_download():
    if exists('ex-ccp.tar.gz'):
        pytest.skip('Data are downloaded.')
    s = 'wget https://osf.io/hzq2x/download -O ex-ccp.tar.gz\n'
    s += 'tar -xzf ex-ccp.tar.gz\n'
    proc = Popen(s, shell=True)
    proc.communicate()


def test_sub01():
    para = ccppara('ex-ccp/ccp.cfg')
    para.rfpath = 'ex-ccp/RFresult'
    para.stalist = 'ex-ccp/sta.lst'
    para.stack_sta_list = ''
    makedata(para)


def test_sub02():
    ccp = CCPProfile('ex-ccp/ccp.cfg')
    # para.stalist = 'ex-ccp/sta.lst'
    ccp.cpara.width = 40
    ccp.initial_profile()
    ccp.stack()
    ccp.save_stack_data(format='dat')
    os.remove(ccp.cpara.stack_sta_list)
    os.remove(ccp.cpara.stackfile)

def test_sub03():
    ccp = CCPProfile('ex-ccp/ccp.cfg')
    ccp.cpara.adaptive = True
    ccp.cpara.stack_sta_list = 'ex-ccp/sta.lst'
    ccp.initial_profile()
    ccp.stack()
    ccp.save_stack_data(format='dat')
    os.remove(ccp.cpara.stackfile)

if __name__ == '__main__':
    # test_download()
    test_sub01()