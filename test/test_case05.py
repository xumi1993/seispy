from seispy.core.depmodel import DepModel
from seispy.seisfwd import SynSeis
from seispy.rfcorrect import RFStation
from seispy.geo import skm2srad
import pytest
import numpy as np


def test_sub01():
    depth = np.array([0, 20.1, 35.1, 100])
    model = DepModel(depth)
    rayp = np.arange(0.04, 0.09, 0.01)
    ss = SynSeis(model, rayp, 0.1, 2400)
    ss.run_fwd()
    rfs = ss.run_deconvolution()
    rfsta = RFStation.read_stream(rfs, rayp, 0)
    rfsta.moveoutcorrect(ref_rayp=0.06, dep_range=np.arange(100))


def test_sub02():
    h = [30, 0]
    vp = [6.5, 8.04]
    vs = [3.6, 4.48]
    rho = [2.67, 3.2]
    rayp = 0.06
    rrayp = skm2srad(rayp)
    dt = 0.1
    depth = np.arange(100)
    basemodel = DepModel.read_layer_model(depth, h, vp, vs, rho=rho)
    basemodel.plot_model(show=False)
    basemodel.tpds(rrayp, rrayp)
    basemodel.tpppds(rrayp, rrayp)
    basemodel.tpspds(rrayp)
    ssfwd = SynSeis(basemodel, rayp, dt, npts=1200)
    ssfwd.run_fwd()


if __name__ == '__main__':
    test_sub01()
