import numpy as np
from seispy import distaz
from pyproj import Geod
from seispy.utils import scalar_instance, array_instance

def sind(deg):
    """ Sine function with degree as input

    :param deg: Degree
    :type deg: float
    :return: Sine value
    :rtype: float
    """
    rad = np.radians(deg)
    return np.sin(rad)


def cosd(deg):
    """ Cosine function with degree as input

    :param deg: Degree
    :type deg: float
    :return: Cosine value
    :rtype: float
    """
    rad = np.radians(deg)
    return np.cos(rad)


def tand(deg):
    """ Tangent function with degree as input

    :param deg: Degree
    :type deg: float
    :return: Tangent value
    :rtype: float
    """
    rad = np.radians(deg)
    return np.tan(rad)


def cotd(deg):
    """ Cotangent function with degree as input

    :param deg: Degree
    :type deg: float
    :return: Cotangent value
    :rtype: float
    """
    rad = np.radians(deg)
    return np.cos(rad) / np.sin(rad)


def asind(x):
    """ Inverse sine function with degree as output

    :param x: Sine value
    :type x: float
    :return: Degree
    :rtype: float
    """
    rad = np.arcsin(x)
    return np.degrees(rad)


def acosd(x):
    """ Inverse cosine function with degree as output

    :param x: Cosine value
    :type x: float
    :return: Degree 
    :rtype: float
    """
    rad = np.arccos(x)
    return np.degrees(rad)


def atand(x):
    """ Inverse tangent function with degree as output

    :param x: Tangent value
    :type x: float
    :return: Degree
    :rtype: float
    """
    rad = np.arctan(x)
    return np.degrees(rad)


def km2deg(km):
    """ Convert km to degree

    :param km: Distance in km
    :type km: float
    :return: Distance in degree
    :rtype: float
    """
    radius = 6371
    circum = 2*np.pi*radius
    conv = circum / 360
    deg = km / conv
    return deg


def deg2km(deg):
    """ Convert degree to km

    :param deg: Distance in degree
    :type deg: float
    :return: Distance in km
    :rtype: float
    """
    radius = 6371
    circum = 2*np.pi*radius
    conv = circum / 360
    km = deg * conv
    return km


def rad2deg(rad):
    """ Convert radian to degree

    :param rad: Radian
    :type rad: float
    :return: Degree
    :rtype: float
    """
    deg = rad*(360/(2*np.pi))
    return deg


def skm2sdeg(skm):
    """ Convert s/km to s/degree

    :param skm: s/km
    :type skm: float
    :return: s/degree
    :rtype: float
    """
    sdeg = skm * deg2km(1)
    return sdeg


def sdeg2skm(sdeg):
    """ Convert s/degree to s/km

    :param sdeg: s/degree
    :type sdeg: float
    :return: s/km
    :rtype: float
    """
    skm = sdeg / deg2km(1)
    return skm


def srad2skm(srad):
    """ Convert s/rad to s/km

    :param srad: s/rad
    :type srad: float
    :return: s/km
    :rtype: float
    """
    sdeg = srad * ((2*np.pi)/360)
    return sdeg / deg2km(1)


def skm2srad(skm):
    """ Convert s/km to s/rad

    :param skm: s/km
    :type skm: float
    :return: s/rad
    :rtype: float
    """
    sdeg = skm * deg2km(1)
    return rad2deg(sdeg)


def rot3D(bazi, inc):
    """
    Rotate components from ZRT to LQT
    
    .. code-block:: python

    M = [[cos(inc)     -sin(inc)*sin(bazi)    -sin(inc)*cos(bazi)],
         [sin(inc)      cos(inc)*sin(bazi)     cos(inc)*cos(bazi)],
         [0              -cos(bazi)             sin(bazi)]]
    
    :param bazi: back-azimuth of station-event pair
    :param inc: Incidence angle
    :return: Coefficient matrix m
    :rtype: np.ndarray
    """

    if scalar_instance(inc):
        value31 = 0
    elif array_instance(inc):
        value31 = np.repeat(0, len(inc))
    else:
        raise TypeError('Input args sould be in \'float\', \'int\', or \'numpy.ndarray\'')

    inc = inc / 180 * np.pi
    bazi = bazi / 180 * np.pi

    m = np.array([[np.cos(inc), -np.sin(inc)*np.sin(bazi), -np.sin(inc)*np.cos(bazi)],
                  [np.sin(inc), np.cos(inc)*np.sin(bazi), np.cos(inc)*np.cos(bazi)],
                  [value31, -np.cos(bazi), np.sin(bazi)]])
    return m


def rotateSeisZENtoLQT(z, e, n, bazi, inc):
    """ Rotate ZEN to LQT

    :param z: Vertical component
    :type z: np.ndarray
    :param e: East component
    :type e: np.ndarray
    :param n: North component
    :type n: np.ndarray
    :param bazi: Back-azimuth
    :type bazi: float
    :param inc: Incidence angle
    :type inc: float
    :return: L, Q and T components
    :rtype: tuple
    """
    m = rot3D(bazi, inc)
    zen = np.array([z, e, n])
    lqt = np.dot(m, zen)
    return lqt[0, :], lqt[1, :], lqt[2, :]


def spherical2cartesian(lon, lat, dep):
    """ Convert spherical coordinates to cartesian coordinates

    :param lon: Longitude
    :type lon: float
    :param lat: Latitude
    :type lat: float
    :param dep: Depth
    :type dep: float
    :return: Cartesian coordinates
    :rtype: tuple
    """
    cola = 90. - lat
    r = 6371 - dep
    x = r * sind(cola) * cosd(lon)
    y = r * sind(cola) * sind(lon)
    z = r * cosd(cola)
    return x, y, z


def rotateSeisENtoTR(e, n, baz):
    """ Rotate EN to TR

    :param e: East component
    :type e: np.ndarray
    :param n: North component
    :type n: np.ndarray
    :param baz: Back-azimuth
    :type baz: float
    :return: T and R components
    :rtype: tuple
    """
    angle = np.mod(baz+180, 360)
    r = n*cosd(angle) + e*sind(angle)
    t = e*cosd(angle) - n*sind(angle)
    return t, r


def rssq(x):
    """ Root sum square
    :param x: Input array
    :type x: np.ndarray
    :return: Root sum square
    :rtype: float
    """
    return np.sqrt(np.sum(x**2)/len(x))


def snr(x, y):
    """ Signal-to-noise ratio
    
    :param x: Signal
    :type x: np.ndarray
    :param y: Noise
    :type y: np.ndarray
    :return: SNR
    :rtype: float
    """
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

    def init_lalo(lat0, lon0, npts):
        if hasattr(lat0, "__iter__") and hasattr(lon0, "__iter__"):
            if len(lat0) != len(lon0):
                raise ValueError('lat0 and lon0 must be in the same length')
            elif len(lat0) != npts:
                raise ValueError('initial points must be in the same length as azimuths')
            else:
                lat1 = lat0
                lon1 = lon0
        elif scalar_instance(lat0) and scalar_instance(lon0):
            lat1 = np.ones(npts) * lat0
            lon1 = np.ones(npts) * lon0
        else:
            raise ValueError('lat0 and lon0 must be in the same length')
        return lat1, lon1

    if hasattr(azimuth, "__iter__") and hasattr(gcarc_dist, "__iter__"):
        if len(azimuth) == len(gcarc_dist):
            npts = len(azimuth)
            lat0, lon0 = init_lalo(lat0, lon0, npts)
        else:
            raise ValueError('azimuth and gcarc_dist must be in the same length')
    elif scalar_instance(azimuth) and scalar_instance(gcarc_dist):
        if hasattr(lat0, "__iter__") and hasattr(lon0, "__iter__"):
            if len(lat0) != len(lon0):
                raise ValueError('lat0 and lon0 must be in the same length')
            else:
                azimuth = np.ones_like(lat0)*azimuth
                gcarc_dist = np.ones_like(lat0)*gcarc_dist
        elif scalar_instance(lat0) and scalar_instance(lon0):
            pass
        else:
            raise ValueError('lat0 and lon0 must be in the same length')            
    elif scalar_instance(azimuth) and hasattr(gcarc_dist, "__iter__"):
        npts = len(gcarc_dist)
        azimuth = np.ones(npts)*azimuth
        lat0, lon0 = init_lalo(lat0, lon0, npts)
    elif scalar_instance(gcarc_dist) and hasattr(azimuth, "__iter__"):
        npts = len(azimuth)
        gcarc_dist = np.ones(lat0, lon0, npts)*gcarc_dist
        lat0, lon0 = init_lalo(lat0, lon0, npts)
    g = Geod(ellps=ellps)
    lon, lat, _ = g.fwd(lon0, lat0, azimuth, deg2km(gcarc_dist)*1000)
    return lat, lon


def geoproject(lat_p, lon_p, lat1, lon1, lat2, lon2):
    """ Project a point to a line

    :param lat_p: Latitude of the point
    :type lat_p: float
    :param lon_p: Longitude of the point
    :type lon_p: float
    :param lat1: Latitude of the first point of the line
    :type lat1: float
    :param lon1: Longitude of the first point of the line
    :type lon1: float
    :param lat2: Latitude of the second point of the line
    :type lat2: float
    :param lon2: Longitude of the second point of the line
    :type lon2: float
    :return: Latitude and longitude of the projected point
    :rtype: tuple
    """
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