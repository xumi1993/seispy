import numpy as np
from numpy import pi, mod
from seispy import distaz
from scipy import interpolate
import math

def sind(deg):
    rad = math.radians(deg)
    return math.sin(rad)

def cosd(deg):
    rad = math.radians(deg)
    return math.cos(rad)

def tand(deg):
    rad = math.radians(deg)
    return math.tan(rad)

def cotd(deg):
    rad = math.radians(deg)
    return math.cos(rad) / math.sin(rad)

def asind(x):
    rad = math.asin(x)
    return math.degrees(rad)

def acosd(x):
    rad = math.acos(x)
    return math.degrees(rad)

def atand(x):
    rad = math.atan(x)
    return math.degrees(rad)

def km2deg(km):
    radius = 6371
    circum = 2*pi*radius
    conv = circum / 360
    deg = km / conv
    return deg

def deg2km(deg):
    radius = 6371
    circum = 2*pi*radius
    conv = circum / 360
    km = deg * conv
    return km

def rad2deg(rad):
    deg = rad*(360/(2*pi))
    return deg

def skm2sdeg(skm):
    sdeg = skm * deg2km(1)
    return sdeg

def sdeg2skm(sdeg):
    skm = sdeg / deg2km(1)
    return skm

def srad2skm(srad):
    sdeg = srad * ((2*pi)/360)
    return sdeg / deg2km(1)


def rot3D(bazi, inc):
    """
    :param bazi:
    :param inc:
    :return:
    M = [cos(inc)     -sin(inc)*sin(bazi)    -sin(inc)*cos(bazi);
        sin(inc)      cos(inc)*sin(bazi)     cos(inc)*cos(bazi);
        0              -cos(bazi)             sin(bazi)];
    """

    if isinstance(inc, float) or isinstance(inc, int):
        value31 = 0
    elif isinstance(inc, np.ndarray):
        value31 = np.repeat(0, len(inc))
    else:
        raise TypeError('Input args sould be in \'float\', \'int\', or \'numpy.ndarray\'')

    inc = inc / 180 * pi
    bazi = bazi / 180 * pi

    M = np.array([[np.cos(inc), -np.sin(inc)*np.sin(bazi), -np.sin(inc)*np.cos(bazi)],
                  [np.sin(inc), np.cos(inc)*np.sin(bazi), np.cos(inc)*np.cos(bazi)],
                  [value31, -np.cos(bazi), np.sin(bazi)]])
    return M


def rotateSeisZENtoLQT(Z, E, N, bazi, inc):
    M = rot3D(bazi, inc)
    ZEN = np.array([Z, E, N])
    LQT = np.dot(M, ZEN)
    return LQT[0, :], LQT[1, :], LQT[2, :]


def rotateSeisENtoTR(E, N, BAZ):
    angle = mod(BAZ+180, 360)
    R = N*cosd(angle) + E*sind(angle)
    T = E*cosd(angle) - N*sind(angle)
    return T, R

def rssq(x):
    return np.sqrt(np.sum(x**2)/len(x))

def snr(x, y):
    spow = rssq(x)**2
    npow = rssq(y)**2
    return 10 * np.log10(spow / npow)

def latlon_from(lat1,lon1,azimuth,gcarc_dist):
    lat2 = asind ((sind (lat1) * cosd (gcarc_dist)) + (cosd (lat1) * sind (gcarc_dist) * cosd (azimuth)))
    if ( cosd (gcarc_dist) >= (cosd (90 - lat1) * cosd (90 - lat2))):
        lon2 = lon1 + asind (sind (gcarc_dist) * sind (azimuth) / cosd (lat2))
    else:
        lon2 = lon1 + asind (sind (gcarc_dist) * sind (azimuth) / cosd (lat2)) + 180
    return lat2, lon2

def geoproject(lat_p, lon_p, lat1, lon1, lat2, lon2):
    azi = distaz(lat1, lon1, lat2, lon2).baz
    dis_center = distaz(lat1, lon1, lat_p, lon_p).delta
    azi_center = distaz(lat1, lon1, lat_p, lon_p).baz
    dis_along = atand(tand(dis_center))*cosd(azi-azi_center)
    (lat, lon) = latlon_from(lat1, lon1, azi, dis_along)
    return lat, lon

def extrema(x, opt='max'):
    if opt == 'max':
        idx = np.intersect1d(np.where(np.diff(x) > 0)[0]+1, np.where(np.diff(x) < 0)[0])
    elif opt == 'min':
        idx = np.intersect1d(np.where(np.diff(x) < 0)[0]+1, np.where(np.diff(x) > 0)[0])
    else:
        raise ImportError('Wrong Options!!!')
    return idx

def slantstack(seis, timeaxis, rayp_range, tau_range, ref_dis, dis):
    EV_num = seis.shape()[1]
    tmp = np.zeros([EV_num, tau_range.shape()[0]])
    amp = np.zeros([rayp_range.shape()[0], tau_range.shape()[0]])
    for j in range(rayp_range.shape()[0]):
        for i in range(EV_num):
            seis[:,i] = seis[:,i] / np.max(np.abs(seis[:,i]))
            tps = tau_range - rayp_range[j] * (dis[i] - ref_dis)
            tmp[i,:] = interpolate.interp1d(timeaxis, seis[:,i].T)(tps)
        amp[j,:] = np.mean(tmp)
    return amp

