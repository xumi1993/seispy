import numpy as np
import argparse
from seispy.rfcorrect import SACStation
from seispy.ccp3d import CCP3D
from seispy.ccpprofile import CCPProfile


def rfani():
    parser = argparse.ArgumentParser(description="Estimate crustal anisotropy with a Joint inversion method. See Liu and Niu (2012) in detail.")
    parser.add_argument('rfpath', type=str, help="Path to PRFs")
    parser.add_argument('-t', help="Time window cut from tb to te", metavar='tb/te', required=True)
    parser.add_argument('-c', help="List file in text format for saving results, defaults to rfani.dat",
                        default="rfani.dat", metavar="list_file_name")
    parser.add_argument('-l', help="Half length of time window when cut out Pms phases",
                        default=3, metavar="half_time_length", type=float)
    parser.add_argument('-o', dest='outpath', help="Directory to the image, defaults to current directory.", default='./')
    parser.add_argument('-p', help="If plot RFs stacked by back-azimuth, defaults to \'False\'",
                        dest="isplot", action='store_true', default=False)
    parser.add_argument('-w', help="Weights of 3 anisotropic methods (order by R cosine energy, R cross-correlation and T energy), defaults to 0.4/0.4/0.2",
                        dest='weight', default='0.4/0.4/0.2', metavar='w1/w2/w3')
    arg = parser.parse_args()
    weights = np.array(arg.weight.split('/')).astype(float)
    timewin = np.array(arg.t.split('/')).astype(float)
    rfsta = SACStation(arg.rfpath)
    bf, bt = rfsta.jointani(timewin[0], timewin[1], tlen=arg.l, weight=weights)
    with open(arg.c, 'a+') as fid:
        for f, t in zip(bf, bt):
            fid.write('{}\t{:.3f}\t{:.3f}\t{:.2f}\t{:.2f}'.format(rfsta.staname, rfsta.stla, rfsta.stlo, f, t))
    if arg.isplot:
        rfsta.ani.plot_stack_baz(outpath=arg.outpath)
    rfsta.ani.plot_polar(outpath=arg.outpath)


def ccp3d():
    parser = argparse.ArgumentParser(description="3-D CCP stacking with spaced grid bins.")
    parser.add_argument('cfg_file', type=str, help='Path to CCP configure file')
    parser.add_argument('-s', help='Range for searching depth of D410 and D660, The results would be saved to \'peakfile\' in cfg_file',
                        metavar='d410min/d410max/d660min/d660max', default=None)
    arg = parser.parse_args()
    ccp = CCP3D(arg.cfg_file)
    ccp.initial_grid()
    ccp.stack()
    ccp.save_stack_data(ccp.cpara.stackfile)
    if arg.s:
        search_range = np.array(arg.s.split('/')).astype(float)
        ccp.search_good_410_660(*search_range)
        ccp.save_good_410_660(ccp.cpara.peakfile)


def ccp_profile():
    parser = argparse.ArgumentParser(description="Stack PRFS along a profile")
    parser.add_argument('cfg_file', type=str, help='Path to CCP configure file')
    parser.add_argument('-t', help='Output as a text file', dest='isdat', action='store_true')
    arg = parser.parse_args()
    if arg.isdat:
        typ = 'dat'
    else:
        typ = 'npz'
    ccp = CCPProfile(arg.cfg_file)
    ccp.initial_profile()
    ccp.stack()
    ccp.save_stack_data(format=typ)
