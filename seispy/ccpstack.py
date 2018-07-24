import numpy as np
from scipy.io import loadmat, savemat
from seispy.geo import km2deg, latlon_from, geoproject
from seispy.ccppara import ccppara
from seispy import distaz
from seispy.rfcorrect import DepModel
from os.path import exists


def gen_profile(lineloca, slid_val):
    da = distaz(lineloca[0, 0], lineloca[0, 1], lineloca[1, 0], lineloca[1, 1])
    profile_range = np.arange(0, da.degreesToKilometers(), slid_val)
    profile_lat, profile_lon = latlon_from(lineloca[0, 0], lineloca[0, 1], da.baz, profile_range)
    return da, profile_range, profile_lat, profile_lon


def read_stations(rfdepth):
    global sta_num
    staname = []
    stla = np.zeros(sta_num)
    stlo = np.zeros(sta_num)
    for i in range(sta_num):
        staname.append(rfdepth[0, i]['Station'][0])
        stla[i] = rfdepth[0, i]['stalat'][0, 0]
        stlo[i] = rfdepth[0, i]['stalon'][0, 0]

    return np.array(staname), stla, stlo


def search_stations(cpara, da, staname, stla, stlo):
    global sta_num
    projstla, projstlo = geoproject(stla, stlo, cpara.line[0, 0], cpara.line[0, 1], cpara.line[1, 0], cpara.line[1, 1])
    sta_dis = distaz(projstla, projstlo, stla, stlo).delta
    baz_sta_start = distaz(cpara.line[0, 0], cpara.line[0, 1], stla, stlo).baz
    baz_sta_end = distaz(cpara.line[1, 0], cpara.line[1, 1], stla, stlo).baz
    start_range = [np.mod(da.az+90, 360), np.mod(da.az-90, 360)]
    end_range = [np.mod(da.baz-90, 360), np.mod(da.baz+90, 360)]
    idx = np.where(((baz_sta_start > start_range[0]) & (baz_sta_start < start_range[1])) & (
                (baz_sta_end < end_range[0]) | (baz_sta_end > end_range[1])) & (sta_dis < cpara.width))[0]
    if idx.size == 0:
        raise ValueError('No station beside this profile')
    if cpara.stack_sta_list != '':
        with open(cpara.stack_sta_list, 'w+') as f:
            for i in idx:
                f.write('{} {:3f} {:3f}\n'.format(staname[i], stla[i], stlo[i]))
    return idx, projstla, projstlo


def gen_bin_radius(cpara):
    if cpara.bin_radius is None:
        velmod = DepModel(cpara.stack_range, cpara.velmod)
        bin_radius = km2deg(np.sqrt(0.5*cpara.domperiod * velmod.vs * cpara.stack_range))
    else:
        bin_radius = np.zeros(cpara.stack_range.shape[0]) * km2deg(cpara.bin_radius)
    return bin_radius


def project(idx, cpara):
    global rfdepth
    dtype = {'names': ('moveout_correct', 'projlat', 'projlon'), 'formats': ('O', 'O', 'O')}
    proj = np.zeros((len(idx), ), dtype=dtype)
    m = 0
    for i in idx:
        ev_eum = rfdepth[0, i]['Piercelat'].shape[1]
        projlat = np.zeros_like(rfdepth[0, i]['Piercelat'])
        projlon = np.zeros_like(rfdepth[0, i]['Piercelon'])
        for j in range(ev_eum):
            projlat[:, j], projlon[:, j] = geoproject(rfdepth[0, i]['Piercelat'][:, j], rfdepth[0, i]['Piercelon'][:, j],
                       cpara.line[0, 0], cpara.line[0, 1], cpara.line[1, 0], cpara.line[1, 1])
        proj[m]['moveout_correct'] = rfdepth[0, i]['moveout_correct']
        proj[m]['projlat'] = projlat
        proj[m]['projlon'] = projlon
        m += 1
    return proj


def find_falling(proj, cpara, profile_range, profile_lat, profile_lon, bin_radius):
    stack_data = np.zeros((profile_range.shape[0], 5), dtype=np.object)
    stack_rf = np.zeros([cpara.stack_range.shape[0], profile_range.shape[0]])
    event_count = np.zeros([cpara.stack_range.shape[0], profile_range.shape[0]])
    dom = int(cpara.stack_val/cpara.dep_val)
    for i in range(profile_range.shape[0]):
        print('calculate the RF stacks at the distance of {}km along the profile-------'.format(profile_range[i]))
        for j in range(cpara.stack_range.shape[0]):
            for k in range(proj.shape[0]):
                col = np.where(distaz(profile_lat[i], profile_lon[i], proj[k]['projlat'][int(dom*cpara.stack_range[j]), :],
                                      proj[k]['projlon'][int(dom*cpara.stack_range[j]), :]).delta < bin_radius[j])[0]
                if col.size != 0:
                    stack_rf[j, i] += np.sum(proj[k]['moveout_correct'][int(dom*cpara.stack_range[j]), col])
                    event_count[j, i] += col.shape[0]
            if event_count[j, i] > 0:
                stack_rf[j, i] /= event_count[j, i]
        stack_data[i, 0] = profile_lat[i]
        stack_data[i, 1] = profile_lon[i]
        stack_data[i, 2] = profile_range[i]
        stack_data[i, 3] = stack_rf[:, i]
        stack_data[i, 4] = event_count[:, i]
    return stack_data


def line_stack(cfg_file):
    global sta_num, rfdepth
    cpara = ccppara(cfg_file)
    rfdepth = loadmat(cpara.depthdat)['RFdepth']
    sta_num = rfdepth.shape[1]
    da, profile_range, profile_lat, profile_lon = gen_profile(cpara.line, cpara.slid_val)
    staname, stla, stlo = read_stations(rfdepth)
    sta_idx, projstla, projstlo = search_stations(cpara, da, staname, stla, stlo)
    bin_radius = gen_bin_radius(cpara)
    proj = project(sta_idx, cpara)

    stack_data = find_falling(proj, cpara, profile_range, profile_lat, profile_lon, bin_radius)
    savemat(cpara.stackfile, {'Stack_data': stack_data}, oned_as='column')


if __name__ == '__main__':
    line_stack('/Users/xumj/Researches/NETibet/Ordos_process/paraCCP.cfg')
