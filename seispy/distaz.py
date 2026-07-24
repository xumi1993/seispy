import numpy as np


# Geographic latitude is converted to geocentric latitude before calculating
# the spherical arc.  This is the factor used by the reference MATLAB
# implementation (distaz.m).
_GEOCENTRIC_FACTOR = 0.993270
_AZIMUTH_ZERO_TOL = 1.0e-5


def _geocentric_latitude(latitude):
    """Return geocentric latitude in radians."""
    latitude = np.deg2rad(latitude)

    # atan2 is equivalent to atan(factor * tan(latitude)) for valid
    # geographic latitudes, but remains well defined at both poles.
    return np.arctan2(
        _GEOCENTRIC_FACTOR * np.sin(latitude),
        np.cos(latitude),
    )


def _normalize_azimuth(angle):
    """Normalize an angle to [0, 360), mapping round-off near 360 to zero."""
    angle = np.mod(np.rad2deg(angle), 360.0)
    near_zero = (
        (np.abs(angle) < _AZIMUTH_ZERO_TOL)
        | (np.abs(angle - 360.0) < _AZIMUTH_ZERO_TOL)
    )
    return np.where(near_zero, 0.0, angle)


def _scalar_or_array(value):
    """Return a float for scalar input and an ndarray for broadcast input."""
    value = np.asarray(value, dtype=float)
    if value.ndim == 0:
        return value.item()
    return value


class distaz:
    """
    Calculate distance, azimuth and back-azimuth between two surface points.

    ``lat1``/``lon1`` describe point 1 (normally the station), and
    ``lat2``/``lon2`` describe point 2 (normally the event). Inputs may be
    scalars or any mutually broadcastable array-like objects.

    The public angle names retain seispy's historical convention:

    * ``baz`` is the bearing from point 1 to point 2 (MATLAB ``daze``).
    * ``az`` is the bearing from point 2 to point 1 (MATLAB ``dazs``).
    * ``delta`` is the geocentric great-circle arc in degrees (MATLAB ``dd``).

    The calculation follows the supplied MATLAB ``distaz.m`` implementation,
    while using division-free ``atan2`` expressions so that poles and
    arbitrary-dimensional NumPy broadcasting are handled safely.

    Parameters
    ----------
    lat1, lon1, lat2, lon2 : float or array-like
        Geographic coordinates in degrees.

    Notes
    -----
    Equations are from Bullen, sections 10.2, pages 154--155. The original
    routine was written by T. Owens (1991) and subsequently ported through
    Fortran, C, Tcl, Java, and Python. NumPy array support was added to seispy
    by Mijian Xu.

    ObsPy's :func:`obspy.geodetics.locations2degrees` uses geographic
    latitudes on a sphere, while this routine first converts them to
    geocentric latitudes to match ``distaz.m``. ObsPy's
    :func:`obspy.geodetics.gps2dist_azimuth` uses a WGS84 ellipsoid by
    default. Results agree when the same geocentric latitudes and a spherical
    Earth (``f=0``) are used.
    """

    def __init__(self, lat1, lon1, lat2, lon2):
        self.stalat = lat1
        self.stalon = lon1
        self.evtlat = lat2
        self.evtlon = lon2

        lat1_array, lon1_array, lat2_array, lon2_array = np.broadcast_arrays(
            np.asarray(lat1, dtype=float),
            np.asarray(lon1, dtype=float),
            np.asarray(lat2, dtype=float),
            np.asarray(lon2, dtype=float),
        )

        geocentric_lat1 = _geocentric_latitude(lat1_array)
        geocentric_lat2 = _geocentric_latitude(lat2_array)

        sin_lat1 = np.sin(geocentric_lat1)
        cos_lat1 = np.cos(geocentric_lat1)
        sin_lat2 = np.sin(geocentric_lat2)
        cos_lat2 = np.cos(geocentric_lat2)

        # Reducing the longitude difference before converting to radians
        # avoids avoidable precision loss for longitudes outside [-180, 180].
        longitude_difference = lon2_array - lon1_array
        wrapped_longitude_difference = (
            np.remainder(longitude_difference + 180.0, 360.0) - 180.0
        )
        longitude_difference_rad = np.deg2rad(wrapped_longitude_difference)
        sin_dlon = np.sin(longitude_difference_rad)
        cos_dlon = np.cos(longitude_difference_rad)

        # The two components below are also the numerator and denominator of
        # the point-1-to-point-2 bearing. Together with the dot product they
        # form atan2(|cross product|, dot product), a stable equivalent of the
        # MATLAB acos/atan distance calculation.
        east = cos_lat2 * sin_dlon
        north = (
            cos_lat1 * sin_lat2
            - sin_lat1 * cos_lat2 * cos_dlon
        )
        dot_product = (
            sin_lat1 * sin_lat2
            + cos_lat1 * cos_lat2 * cos_dlon
        )
        cross_product_norm = np.hypot(east, north)

        delta = np.rad2deg(np.arctan2(cross_product_norm, dot_product))
        baz = _normalize_azimuth(np.arctan2(east, north))

        # Reverse bearing: point 2 to point 1. This is MATLAB's ``dazs`` and
        # seispy's historical ``az`` attribute.
        reverse_east = -cos_lat1 * sin_dlon
        reverse_north = (
            cos_lat2 * sin_lat1
            - sin_lat2 * cos_lat1 * cos_dlon
        )
        az = _normalize_azimuth(np.arctan2(reverse_east, reverse_north))

        # Longitudes differing by full rotations identify the same point.
        # At either pole, longitude is immaterial. Bearings for coincident
        # points are undefined, so retain distaz.m's established zero value.
        same_latitude = lat1_array == lat2_array
        same_longitude = np.remainder(longitude_difference, 360.0) == 0.0
        same_pole = same_latitude & (np.abs(lat1_array) == 90.0)
        coincident = same_latitude & (same_longitude | same_pole)

        delta = np.where(coincident, 0.0, delta)
        az = np.where(coincident, 0.0, az)
        baz = np.where(coincident, 0.0, baz)

        self.delta = _scalar_or_array(delta)
        self.az = _scalar_or_array(az)
        self.baz = _scalar_or_array(baz)

    def getDelta(self):
        return self.delta

    def getAz(self):
        return self.az

    def getBaz(self):
        return self.baz

    def degreesToKilometers(self):
        return self.delta * 111.19
