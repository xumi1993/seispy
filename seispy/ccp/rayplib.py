"""
use obspy.taup instead of taup in java
"""

import numpy as np

from obspy.taup import TauPyModel

def rayp_input(phase_main:str,
               source_dist_in_degree:str,
               conver_layers:str,
               source_depth_in_km:str,
               out_put_path:str):
    """
    anaylize input informations
    only allow P, SKS, S, ScS
    """
    if phase_main not in ['P', 'S', 'SKS', 'ScS']:
        raise ValueError('not a valid phase')

    try:
        # set default step value for 2 input
        source_dist_in_degree +='/1';source_dist_in_degree = [float(_i) for _i in source_dist_in_degree.split('/')[:4]]
        conver_layers +='/1';conver_layers = [float(_i) for _i in conver_layers.split('/')[:4]]
        source_depth_in_km +='/2';source_depth_in_km = [float(_i) for _i in source_depth_in_km.split('/')]
        if conver_layers[0] < 0:
            conver_layers[0] = 0
        conver_layers = np.arange(conver_layers[0], conver_layers[1], conver_layers[2])
        source_dist_in_degree = np.arange(source_dist_in_degree[0], source_dist_in_degree[1], source_dist_in_degree[2])
        source_depth_in_km = np.arange(source_depth_in_km[0],source_depth_in_km[1], source_depth_in_km[2])
    except:
        raise ValueError('invalid inputs in dist and layers')
    # here we limit layer_max to 30km
    if conver_layers[-1] < 30:
        raise ValueError('theres no need to do such a thing')


    rayp = cal_taup('iasp91', phase_main,
                    conver_layers, source_dist_in_degree, source_depth_in_km)

    np.savez(out_put_path, dis = source_dist_in_degree, dep = source_depth_in_km, layers = conver_layers, rayp = rayp)




def cal_taup(model:str, phase:str,
             layer:np.ndarray, dist:np.ndarray, depth:np.ndarray):
    """
    cal rayp from obspy.taup
    wont cal layers above 11km, change min_con

    phase_list: list of phase names, like ["P30s", "P50s"]
    layer : conversion layers
    dist : numpy.ndarray, dist in deg
    depth : numpy.ndarray, source_depth in km

    """
    min_con = 11
    if phase[-1] == 'P':
        conv = 's'
    else:
        conv = 'p'
    ## init taup and matrix
    rayp_lib = np.zeros((dist.shape[0], depth.shape[0], layer.shape[0]))
    model = TauPyModel(model=model)
    cal_layers_slice = np.where(layer > min_con); head = cal_layers_slice[0][0]
    phase_list = ["{}{}{}".format(phase,_d, conv) for _d in layer[cal_layers_slice]]

    ### main
    for _i, _deg in enumerate(dist):
        for _j, _dep in enumerate(depth):
            print("now at degree {}, dep{}".format(_deg, _dep))
            arrivals=model.get_travel_times(source_depth_in_km=_dep,
                                            distance_in_degree=_deg,
                                            phase_list=phase_list)
            for _k in range(len(arrivals)):
                # 101 for change rayp from s/deg to s/km
                rayp_lib[_i,_j,head+_k] = arrivals[_k].ray_param_sec_degree
    if head > 0:
        for _layer in range(head):
            rayp_lib[:,:,_layer] = rayp_lib[:,:,head]
    return rayp_lib


if __name__ == "__main__":
    phase = "P"
    dist = "30/50/1"
    layers = "0/100/2"
    depth = '0/100/5'
    out_put = "G:/Ps_rayp.npz"
    rayp_input("P", dist, layers, depth, out_put)
