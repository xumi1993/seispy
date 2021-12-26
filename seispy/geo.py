import numpy as np
from numpy import pi, mod
from seispy import distaz

def sind(deg):
    rad = np.radians(deg)
    return np.sin(rad)


def cosd(deg):
    rad = np.radians(deg)
    return np.cos(rad)


def tand(deg):
    rad = np.radians(deg)
    return np.tan(rad)


def cotd(deg):
    rad = np.radians(deg)
    return np.cos(rad) / np.sin(rad)


def asind(x):
    rad = np.arcsin(x)
    return np.degrees(rad)


def acosd(x):
    rad = np.arccos(x)
    return np.degrees(rad)


def atand(x):
    rad = np.arctan(x)
    return np.degrees(rad)


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


def skm2srad(skm):
    sdeg = skm * deg2km(1)
    return rad2deg(sdeg)


def rot3D(bazi, inc):
    """
    Rotate components from ZRT to LQT
    M = [cos(inc)     -sin(inc)*sin(bazi)    -sin(inc)*cos(bazi);
        sin(inc)      cos(inc)*sin(bazi)     cos(inc)*cos(bazi);
        0              -cos(bazi)             sin(bazi)];
    
    :param bazi: back-azimuth of station-event pair
    :param inc: Incidence angle
    :return: Coefficient matrix M
    
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


def spherical2cartesian(lon, lat, dep):
    cola = 90. - lat
    r = 6371 - dep
    x = r * sind(cola) * cosd(lon)
    y = r * sind(cola) * sind(lon)
    z = r * cosd(cola)
    return x, y, z


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
    if npow == 0:
        npow = 0.001
    return 10 * np.log10(spow / npow)


def latlon_from(lat1, lon1, azimuth, gcarc_dist):
    lat2 = asind((sind(lat1) * cosd(gcarc_dist)) + (cosd(lat1) * sind(gcarc_dist) * cosd(azimuth)))
    if isinstance(gcarc_dist, np.ndarray):
        lon2 = np.zeros_like(lat2)
        for n in range(len(gcarc_dist)):
            if cosd(gcarc_dist[n]) >= (cosd(90 - lat1) * cosd(90 - lat2[n])):
                lon2[n] = lon1 + asind(sind(gcarc_dist[n]) * sind(azimuth) / cosd(lat2[n]))
            else:
                lon2[n] = lon1 + asind(sind(gcarc_dist[n]) * sind(azimuth) / cosd(lat2[n])) + 180
    elif isinstance(azimuth, np.ndarray):
        lon2 = np.zeros_like(lat2)
        for n in range(len(azimuth)):
            if cosd(gcarc_dist) >= (cosd(90 - lat1) * cosd(90 - lat2[n])):
                lon2[n] = lon1 + asind(sind(gcarc_dist) * sind(azimuth[n]) / cosd(lat2[n]))
            else:
                lon2[n] = lon1 + asind(sind(gcarc_dist) * sind(azimuth[n]) / cosd(lat2[n])) + 180
    else:
        if (cosd(gcarc_dist) >= (cosd(90 - lat1) * cosd(90 - lat2))):
            lon2 = lon1 + asind(sind(gcarc_dist) * sind(azimuth) / cosd(lat2))
        else:
            lon2 = lon1 + asind(sind(gcarc_dist) * sind(azimuth) / cosd(lat2)) + 180
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
        raise ValueError('opt must be \'max\' or \'min\'')
    return idx


def geo2sph(dep, lat, lon):
    """Convert geographic coordinates to spherical coordinates.

    :param dep: Depth from surface
    :type dep: float or np.ndarray
    :param lat: Latitude
    :type lat: float or np.ndarray
    :param lon: Longitude
    :type lon: float or np.ndarray
    :return: spherical coordinate r, theta, phi
    :rtype: list
    """
    theta = np.radians(90. - lat)
    phi = np.radians(lon)
    r = 6371 - dep
    return r, theta, phi

def sph2geo(r, theta, phi):
    """
    Convert spherical coordinates to geographic coordinates.
    :param float r: radial distance from coordinate system origin
                    {**Units**: km, **Range**: [0, inf)}
    :param float theta: polar angle {**Units**: radians, **Range**: [0,
                        π]}
    :param float phi: azimuthal angle {**Units**: radians, **Range**:
                      [-π, π]}
    :returns: geographic coordinate conversion *(lat, lon, depth)* of
              spherical coordinates
    :rtype: (float, float, float)
    """
    z = 6371 - r
    lat = 90 - rad2deg(theta)
    lon = rad2deg(phi)
    return z, lat, lon