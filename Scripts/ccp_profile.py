import numpy as np
from seispy.geo import distaz, km2deg, deg2km, latlon_from
from seispy.bootstrap import ci
from seispy.distaz import distaz


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


def search_pierce(rfdep, depaxis, bin_loca, profile_range, bin_radius=75):
    bin_radius = km2deg(bin_radius)
    data = []
    for i in range(bin_loca.shape[0]):
        rfbin = {}
        print('{}/{}'.format(i, bin_loca.shape[0]))
        ccp_mean = np.zeros(depaxis.shape[0])
        ccp_count = np.zeros(depaxis.shape[0])
        ccp_ci = np.zeros((depaxis.shape[0], 2))
        for j, dep in zip(range(depaxis.shape[0]), depaxis):
            bin_dep = np.array([])
            for sta in rfdep:
                fall_idx = np.where(distaz(sta['Piercelat'][dep, :], sta['Piercelon'][dep, :], bin_loca[i, 0],
                                           bin_loca[i, 1]).delta < bin_radius)[0]
                bin_dep = np.append(bin_dep, sta['moveout_correct'][dep, fall_idx])
            if bin_dep.shape[0] > 1:
                bin_ci = ci(bin_dep, n_samples=2000)
                bin_mu = np.average(bin_dep)
            else:
                bin_ci = (np.nan, np.nan)
                bin_mu = np.nan
            ccp_count[j] = bin_dep.shape[0]
            ccp_mean[j] = bin_mu
            ccp_ci[j, 0] = bin_ci[0]
            ccp_ci[j, 1] = bin_ci[1]
        rfbin['bin_lat'] = bin_loca[i, 0]
        rfbin['bin_lon'] = bin_loca[i, 1]
        rfbin['profile_dis'] = profile_range[i]
        rfbin['mu'] = ccp_mean
        rfbin['ci'] = ccp_ci
        rfbin['count'] = ccp_count
        data.append(rfbin)
    return data


if __name__ == '__main__':
    lat1 = 37.8
    lon1 = 96
    lat2 = 25.8
    lon2 = 92
    bin_loca, profile_range = init_profile(lat1, lon1, lat2, lon2, 25)
    rfdep = np.load('/Users/xumj/Researches/Tibet_MTZ/RFdepth_1D.npy')
    depaxis = np.arange(300, 800)
    ccp_data = search_pierce(rfdep, depaxis, bin_loca, profile_range)
    np.save('/Users/xumj/Researches/Tibet_MTZ/ccp_B', ccp_data)