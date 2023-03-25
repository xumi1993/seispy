import numpy as np
from numpy import pi, mod
from seispy import distaz
from pyproj import Geod

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


def latlon_from(lat0, lon0, azimuth, gcarc_dist, ellps="WGS84"):
    """
    Determine position with given position of initial point, azimuth and distance

    Accepted numeric scalar or array:
    - :class:`int`
    - :class:`float`
    - :class:`numpy.floating`
    - :class:`numpy.integer`
    - :class:`list`
    - :class:`tuple`
    - :class:`array.array`
    - :class:`numpy.ndarray`
    - :class:`xarray.DataArray`
    - :class:`pandas.Series`

    :param lat0: Latitude of original point
    :type lat0: float or array
    :param lon0: Longitude of original point
    :type lon0: float or array
    :param azimuth: Azimuth(s) in degree
    :type azimuth: float or array
    :param gcarc_dist: Distance(s) between initial and terminus point(s) in degree
    :type gcarc_dist: float or array
    :param ellps: Ellipsoids supported by ``pyproj``, defaults to "WGS84"
    :type ellps: :class:`str`, optional

    Returns
    -------
    scalar or array:
        Latitude(s) of terminus point(s)  
    scalar or array:
        Longitude(s) of terminus point(s)
    """

    if hasattr(azimuth, "__iter__") and hasattr(gcarc_dist, "__iter__"):
        if len(azimuth) == len(gcarc_dist):
            npnt = len(azimuth)
            if hasattr(lat0, "__iter__") and hasattr(lon0, "__iter__"):
                if len(lat0) != len(lon0):
                    raise ValueError('lat0 and lon0 must be in the same length')
                elif len(lat0) != len(azimuth):
                    raise ValueError('initial points must be in the same length as azimuths')
            elif isinstance(lat0, (int, float)) and isinstance(lon0, (int, float)):
                lat0 = np.ones(npnt) * lat0
                lon0 = np.ones(npnt) * lon0
            else:
                raise ValueError('lat0 and lon0 must be in the same length')
        else:
            raise ValueError('azimuth and gcarc_dist must be in the same length')
    elif isinstance(azimuth, (int, float)) and isinstance(gcarc_dist, (int, float)):
        if hasattr(lat0, "__iter__") and hasattr(lon0, "__iter__"):
            if len(lat0) != len(lon0):
                raise ValueError('lat0 and lon0 must be in the same length')
            else:
                azimuth = np.ones(lat0)*azimuth
                gcarc_dist = np.ones(lat0)*gcarc_dist
        elif isinstance(lat0, (int, float)) and isinstance(lon0, (int, float)):
            pass
        else:
            raise ValueError('lat0 and lon0 must be in the same length')            

    g = Geod(ellps=ellps)
    lon, lat, _ = g.fwd(lon0, lat0, azimuth, deg2km(gcarc_dist)*1000)
    return lat, lon


def geoproject(lat_p, lon_p, lat1, lon1, lat2, lon2):
    azi = distaz(lat1, lon1, lat2, lon2).baz
    dis_center = distaz(lat1, lon1, lat_p, lon_p).delta
    azi_center = distaz(lat1, lon1, lat_p, lon_p).baz
    dis_along = atand(tand(dis_center))*cosd(azi-azi_center)
    lat, lon = latlon_from(lat1, lon1, azi, dis_along)
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