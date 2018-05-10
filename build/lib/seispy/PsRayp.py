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
    arr_num = len(arrs)
    print(arr_num, len(ph_list))
    rayp_list = np.zeros([arr_num, 2])
    for i in range(arr_num):
        rayp_list[i][0] = srad2skm(arrs[i].ray_param)
        rayp_list[i][1] = int(arrs[i].name.strip('Ps'))
    rayp_list.sort(axis=0)
    return(rayp_list)


if __name__ == '__main__':
    layers = np.arange(1,700)
    print(PsRayp(layers, 50, 10))
