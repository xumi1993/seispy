import numpy as np
from seispy.geo import distaz, km2deg
from seispy.bootstrap import ci
from scipy.io import savemat, loadmat
import warnings


warnings.filterwarnings("ignore")


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


def boot_stack(ccp_data, bin_loca, depaxis):
    stack_data = []
    for i in range(ccp_data.shape[0]):
        print('{}/{}'.format(i+1, bin_loca.shape[0]))
        boot_stack = {}
        bin_mu = np.zeros(depaxis.shape[0])
        bin_ci = np.zeros([depaxis.shape[0], 2])
        bin_count = np.zeros(depaxis.shape[0])
        for j in range(depaxis.shape[0]):
            bin_count[j] = ccp_data[i, j].shape[1]
            if ccp_data[i, j].shape[1] > 1:
                cci = ci(ccp_data[i, j], n_samples=2000)
                bin_ci[j, 0] = cci[0]
                bin_ci[j, 1] = cci[1]
                bin_mu[j] = np.average(ccp_data[i, j])
            else:
                bin_ci[j, 0] = np.nan
                bin_ci[j, 1] = np.nan
                bin_mu[j] = np.nan
        boot_stack['bin_lat'] = bin_loca[i, 0]
        boot_stack['bin_lon'] = bin_loca[i, 1]
        boot_stack['mu'] = bin_mu
        boot_stack['ci'] = bin_ci
        boot_stack['count'] = bin_count
        stack_data.append(boot_stack)
    return stack_data


if __name__ == '__main__':
    lat1 = 25
    lon1 = 81
    lat2 = 39
    lon2 = 101
    bin_loca = initgrid(lat1, lon1, lat2, lon2, 0.5)
    # rfdep = np.load('/Volumes/xumj3/TibetRF/RFdepth_1D.npy')
    depaxis = np.arange(300, 800)
    # data = search_pierce(rfdep, depaxis, bin_loca)
    ccp_data = loadmat('/Users/xumj/Researches/Tibet_MTZ/ccp_data.mat')['ccp_data']
    stack_data = boot_stack(ccp_data, bin_loca, depaxis)
    np.save('/Users/xumj/Researches/Tibet_MTZ/ccp_results/stack_data', stack_data)