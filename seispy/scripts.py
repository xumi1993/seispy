import numpy as np
import argparse
from seispy.rfcorrect import SACStation


def rfani():
    parser = argparse.ArgumentParser(description="Estimate crustal anisotropy with a Joint inversion method. See Liu and Niu (2012) in detail.")
    parser.add_argument('rfpath', type=str, help="Path to PRFs")
    parser.add_argument('-c', help="If the option was specified directly output results without commentary", action='store_true', default=False)
    parser.add_argument('-o', dest='outpath', help="Directory to the image.")
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