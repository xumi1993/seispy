import numpy as np
import argparse
from seispy.rfcorrect import SACStation
from seispy.ccp3d import CCP3D
from seispy.ccpprofile import CCPProfile


def rfani():
    parser = argparse.ArgumentParser(description="Estimate crustal anisotropy with a Joint inversion method. See Liu and Niu (2012) in detail.")
    parser.add_argument('rfpath', type=str, help="Path to PRFs")
    parser.add_argument('-c', help="If the option was specified, directly output results without commentary", action='store_true', default=False)
    parser.add_argument('-o', dest='outpath', help="Directory to the image.", default='./')
    parser.add_argument('-p', help="Plot energy of anisotropy and show via Matplotlib",
                        dest="isplot", action='store_true', default=False)
    parser.add_argument('-t', help="Time window cut from tb to te", metavar='tb/te', required=True)
    parser.add_argument('-w', help="Weights of 3 anisotropic methods (order by R cosine energy, R cross-correlation and T energy)", 
                        dest='weight', default='0.4/0.4/0.2', metavar='w1/w2/w3')
    arg = parser.parse_args()
    weights = np.array(arg.weight.split('/')).astype(float)
    timewin = np.array(arg.t.split('/')).astype(float)
    rfsta = SACStation(arg.rfpath)
    bf, bt = rfsta.jointani(timewin[0], timewin[1], weight=weights)
    if not arg.c:
        header = 'FVD\tdt\n-----------'
        print(header)
    for f, t in zip(bf, bt):
        print('{}\t{}'.format(f, t))
    rfsta.ani.plot_polar(show=arg.isplot, outpath=arg.outpath)


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
    # parser.add_argument('-d', help='Path to 3d vel model in npz file', dest='vel3dpath', type=str, default='')
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

    