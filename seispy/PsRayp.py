import obspy
from obspy.taup import TauPyModel
import numpy as np
from seispy.geo import srad2skm

def makepheaselist(layers):
    ph_list = ['P'+str(depth)+'s' for depth in layers]
    return ph_list

def PsRayp(layers, dist, dep):
    model = TauPyModel(model="iasp91")
    ph_list = makepheaselist(layers)
    arrs = model.get_ray_paths(dist, dep, ph_list)
    return [(srad2skm(arr.ray_param), arr.name) for arr in arrs]

if __name__ == '__main__':
    layers = np.arange(1,700)
    print(PsRayp(layers, 50, 10))
