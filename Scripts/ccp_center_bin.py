import numpy as np
from seispy.geo import distaz, km2deg


def initgrid(lat1, lon1, lat2, lon2, val):
    LonRange = np.arange(lon1, lon2 + val, val)
    LatRange = np.arange(lat1, lat2 + val, val)
    bin_loca = np.zeros((LonRange.shape[0]*LatRange.shape[0], 2))
    i = 0
    for lat in LatRange:
        for lon in LonRange:
            bin_loca[i, 0] = lat
            bin_loca[i, 1] = lon
            i += 1
    return bin_loca


def search_pierce(rfdep, depaxis, bin_loca, bin_radius=100):
    bin_radius = km2deg(bin_radius)
    data = np.zeros((bin_loca.shape[0], depaxis.shape[0]), dtype='O')
    for i in range(bin_loca.shape[0]):
        print('{}/{}'.format(i, bin_loca.shape[0]))
        for j, dep in zip(range(depaxis.shape[0]), depaxis):
            bin_dep = np.array([])
            for sta in rfdep:
                fall_idx = np.where(distaz(sta['Piercelat'][dep, :], sta['Piercelon'][dep, :], bin_loca[i, 0],
                                           bin_loca[i, 1]).delta < bin_radius)[0]
                bin_dep = np.append(bin_dep, sta['moveout_correct'][dep, fall_idx])
            data[i, j] = bin_dep
    return data


if __name__ == '__main__':
    lat1 = 25
    lon1 = 81
    lat2 = 39
    lon2 = 101
    bin_loca = initgrid(lat1, lon1, lat2, lon2, 0.5)
    rfdep = np.load('/Volumes/xumj3/TibetRF/RFdepth_1D.npy')
    depaxis = np.arange(300, 800)
    data = search_pierce(rfdep, depaxis, bin_loca)
    np.savez('/Volumes/xumj3/TibetRF/ccp_data', ccp_data=data)
