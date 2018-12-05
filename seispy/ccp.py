import numpy as np
from seispy.geo import km2deg, deg2km, latlon_from, geoproject
from seispy.distaz import distaz
from seispy.ccppara import ccppara
from seispy.setuplog import setuplog
from seispy.rfcorrect import DepModel
from seispy.bootstrap import ci
from scipy.interpolate import interp1d
from scipy.io import loadmat
import argparse
import sys
from os.path import dirname


def init_profile(lat1, lon1, lat2, lon2, val):
    azi = distaz(lat1, lon1, lat2, lon2).baz
    dis = distaz(lat1, lon1, lat2, lon2).delta
    profile_range = np.arange(0, deg2km(dis), val)
    lat_loca, lon_loca = latlon_from(lat1, lon1, azi, km2deg(profile_range))
    bin_loca = np.zeros([lat_loca.shape[0], 2])
    for i in range(lat_loca.shape[0]):
        bin_loca[i, 0] = lat_loca[i]
        bin_loca[i, 1] = lon_loca[i]
    return bin_loca, profile_range


def get_sta(rfdep, stalist, line_loca, dep_axis, log):
    new_rfdep = []
    sta_name_all = [st['Station'][0, 0][0] for st in rfdep]
    log.CCPlog.info('Use stacking stations in {}'.format(stalist))
    sta_name = list(np.loadtxt(stalist, dtype=np.dtype('U10'), usecols=(0,), ndmin=1, unpack=True))
    # with open(stalist) as f:
    #     sta_name = [line.split()[0] for line in f.readlines()]
    for sta in sta_name:
        try:
            idx = sta_name_all.index(sta)
        except ValueError:
            log.CCPlog.warning('{} does not in RFdepth structure'.format(sta))
            continue
        projlat = np.zeros_like(rfdep[idx]['Piercelat'][0, 0])
        projlon = np.zeros_like(rfdep[idx]['Piercelon'][0, 0])
        for i, dep in enumerate(dep_axis):
            projlat[:, i], projlon[:, i] = geoproject(rfdep[idx]['Piercelat'][0, 0][:, i],
                                                          rfdep[idx]['Piercelon'][0, 0][:, i],
                                                          line_loca[0, 0], line_loca[0, 1],
                                                          line_loca[1, 0], line_loca[1, 1])
        this_dtype = np.dtype(rfdep[idx].dtype.descr+[('projlat', 'f8', projlat.shape), ('projlon', 'f8', projlon.shape)])
        this_rfdep = np.empty(rfdep[idx].shape, this_dtype)
        for field in rfdep[idx].dtype.descr:
            this_rfdep[field[0]] = rfdep[idx][field[0]]
        this_rfdep['projlat'] = projlat
        this_rfdep['projlon'] = projlon
        new_rfdep.append(this_rfdep)
    return new_rfdep


def search_pierce(rfdep, bin_loca, profile_range, stack_range, dep_axis, log, bin_radius=None, isci=False, domperiod=5):
    dep_idx = np.arange(dep_axis.shape[0])
    stack_idx = np.int16(np.round(interp1d(dep_axis, dep_idx)(stack_range)))
    data = []
    depmod = DepModel(dep_axis)
    fzone = km2deg(np.sqrt(0.5*domperiod*depmod.vs*dep_axis))
    for i in range(bin_loca.shape[0]):
        rfbin = {}
        log.CCPlog.info('{}/{} bin from {} at lat: {} lon: {}'.format(i + 1, bin_loca.shape[0], profile_range[i],
                                                                      bin_loca[i, 0], bin_loca[0, 1]))
        ccp_mean = np.zeros(stack_range.shape[0])
        ccp_count = np.zeros_like(ccp_mean)
        if isci:
            ccp_ci = np.zeros((stack_range.shape[0], 2))
        for jj, j in enumerate(stack_idx):
            if bin_radius is None:
                bin_radius = fzone[j]
            bin_dep = np.array([])
            for sta in rfdep:
                fall_idx = np.where(distaz(sta['projlat'][0, 0][:, j], sta['projlon'][0, 0][:, j],
                                           bin_loca[i, 0], bin_loca[i, 1]).delta < bin_radius)[0]
                bin_dep = np.append(bin_dep, sta['moveout_correct'][0, 0][fall_idx, j])
            if bin_dep.shape[0] > 1:
                if isci:
                    ccp_ci[jj] = ci(bin_dep, n_samples=2000)
                ccp_count[jj] = bin_dep.shape[0]
                ccp_mean[jj] = np.average(bin_dep)
            else:
                ccp_mean[jj] = np.nan
        rfbin['bin_lat'] = bin_loca[i, 0]
        rfbin['bin_lon'] = bin_loca[i, 1]
        rfbin['profile_dis'] = profile_range[i]
        rfbin['mu'] = ccp_mean
        if isci:
            rfbin['ci'] = ccp_ci
        rfbin['count'] = ccp_count
        data.append(rfbin)
    return data


def stack(rfdep, cpara, log=setuplog()):
    """
    :param rfdep: RFdepth struct
    :param cpara: parameters for CCP stacking
    :param log: class for seispy.setuplog.setuplog
    :return:
    """
    bin_radius = km2deg(cpara.bin_radius)
    bin_loca, profile_range = init_profile(cpara.line[0, 0], cpara.line[0, 1], cpara.line[1, 0],
                                           cpara.line[1, 1], cpara.slid_val)
    new_rfdep = get_sta(rfdep, cpara.stack_sta_list, cpara.line, cpara.depth_axis, log)
    # del rfdep
    stack_data = search_pierce(new_rfdep, bin_loca, profile_range,
                               cpara.stack_range, cpara.depth_axis, log, bin_radius=bin_radius, isci=False, domperiod=5)
    return stack_data


def ccp_profile():
    parser = argparse.ArgumentParser(description="Stack PRFS along a profile")
    # parser.add_argument('-d', help='Path to 3d vel model in npz file', dest='vel3dpath', type=str, default='')
    parser.add_argument('cfg_file', type=str, help='Path to CCP configure file')
    parser.add_argument('-l', help='Location of the profile <lat1>/<lon1>/<lat2>/<lon2>', dest='line_loca', type=str, default=None)
    arg = parser.parse_args()
    log = setuplog()
    cpara = ccppara(arg.cfg_file)
    if arg.line_loca is not None:
        try:
            line_loca = [float(value) for value in arg.line_loca.split('/')]
        except Exception:
            log.CCPlog.error('Error in \'-l\' option. The format should be <lat1>/<lon1>/<lat2>/<lon2>')
            sys.exit(1)
        cpara.line = np.array([[line_loca[0], line_loca[1]], [line_loca[2], line_loca[3]]])
    rfdep = loadmat(cpara.depthdat)['RFdepth'][0, :]
    try:
        stack_data = stack(rfdep, cpara, log)
    except Exception as e:
        log.CCPlog.error('{}'.format(e))
        sys.exit(1)

    try:
        np.save(cpara.stackfile, stack_data)
    except FileNotFoundError:
        log.RF2depthlog.warning('No such file or directory: {}'.format(dirname(cpara.depthdat)))
        stack_path = input('Enter a exist path:')
        np.save(stack_path, stack_data)


if __name__ == '__main__':
    cpara = ccppara('/Users/xumj/Researches/NETibet/Ordos_Process/paraCCP.cfg')
    rfdep = loadmat(cpara.depthdat)['RFdepth'][0, :]
    stack_data = stack(rfdep, cpara)
    np.save(cpara.stackfile, stack_data)