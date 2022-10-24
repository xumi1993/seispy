from seispy.utils import DepModel
from seispy.seisfwd import SynSeis
from seispy.rfcorrect import RFStation
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


if __name__ == '__main__':
    test_sub01()