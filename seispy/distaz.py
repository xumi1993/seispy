import numpy as np


class distaz:
    """
    DistAZ class

    Calculate the distance, azimuth and back-azimuth between two points on the
    Earth's surface.

    :param lat1: Latitude of point 1
    :type lat1: float
    :param lon1: Longitude of point 1
    :type lon1: float
    :param lat2: Latitude of point 2
    :type lat2: float or np.ndarray
    :param lon2: Longitude of point 2
    :type lon2: float or np.ndarray
    :return: An instance of DistAZ
    :rtype: DistAZ

     Subroutine to calculate the Great Circle Arc distance
        between two sets of geographic coordinates
    
     Equations take from Bullen, pages 154, 155
    
     T. Owens, September 19, 1991
               Sept. 25 -- fixed az and baz calculations
    
    P. Crotwell, Setember 27, 1995
    Converted to c to fix annoying problem of fortran giving wrong
    answers if the input doesn't contain a decimal point.
    
    H. P. Crotwell, September 18, 1997
    Java version for direct use in java programs.
    *
    * C. Groves, May 4, 2004
    * Added enough convenience constructors to choke a horse and made public double
    * values use accessors so we can use this class as an immutable
    H.P. Crotwell, May 31, 2006
    Port to python, thus adding to the great list of languages to which
    distaz has been ported from the origin fortran: C, Tcl, Java and now python
    and I vaguely remember a perl port. Long live distaz! 
    Mijian Xu, Jan 01, 2016: Add compatibility of np.ndarray for input variables.
    Mijian XU, Jun 12, 2025: Update to latest numpy standards.

    """

    def __init__(self, lat1, lon1, lat2, lon2):
        
        self.stalat = lat1
        self.stalon = lon1
        self.evtlat = lat2
        self.evtlon = lon2
        '''
        if (lat1 == lat2) and (lon1 == lon2):
            self.delta = 0.0
            self.az = 0.0
            self.baz = 0.0
            return
        '''

        rad = 2. * np.pi / 360.0
        '''
        scolat and ecolat are the geocentric colatitudes
        as defined by Richter (pg. 318)
        
        Earth Flattening of 1/298.257 take from Bott (pg. 3)
        
        '''
        sph = 1.0 / 298.257

        scolat = np.pi / 2.0 - np.arctan((1. - sph) * (1. - sph) * np.tan(lat1 * rad))
        ecolat = np.pi / 2.0 - np.arctan((1. - sph) * (1. - sph) * np.tan(lat2 * rad))
        slon = lon1 * rad
        elon = lon2 * rad
        """
    
	    a - e are as defined by Bullen (pg. 154, Sec 10.2)
	    These are defined for the pt. 1

        """
        a = np.sin(scolat) * np.cos(slon)
        b = np.sin(scolat) * np.sin(slon)
        c = np.cos(scolat)
        d = np.sin(slon)
        e = -np.cos(slon)
        g = -c * e
        h = c * d
        k = -np.sin(scolat)
        """
	c
	c  aa - ee are the same as a - e, except for pt. 2
	c
        """
        aa = np.sin(ecolat) * np.cos(elon)
        bb = np.sin(ecolat) * np.sin(elon)
        cc = np.cos(ecolat)
        dd = np.sin(elon)
        ee = -np.cos(elon)
        gg = -cc * ee
        hh = cc * dd
        kk = -np.sin(ecolat)
        """
	c
	c  Bullen, Sec 10.2, eqn. 4
	c
        """
        clamped_value = np.clip(a * aa + b * bb + c * cc, -1.0, 1.0)
        delrad = np.arccos(clamped_value)
        self.delta = delrad / rad
        """
	c
	c  Bullen, Sec 10.2, eqn 7 / eqn 8
	c
	c    pt. 1 is unprimed, so this is technically the baz
	c
	c  Calculate baz this way to avoid quadrant problems
	c
        """
        rhs1 = (aa - d) * (aa - d) + (bb - e) * (bb - e) + cc * cc - 2.
        rhs2 = (aa - g) * (aa - g) + (bb - h) * (bb - h) + (cc - k) * (cc - k) - 2.
        dbaz = np.arctan2(rhs1, rhs2)

        # dbaz_idx = np.where(dbaz < 0.0)[0]
        dbaz_idx = np.atleast_1d(dbaz < 0.0).nonzero()[0]
        if len(dbaz_idx) != 0:
            if np.isscalar(dbaz):
                dbaz += 2 * np.pi
            else:
                dbaz[dbaz_idx] += 2 * np.pi

        self.baz = dbaz / rad
        """
	c
	c  Bullen, Sec 10.2, eqn 7 / eqn 8
	c
	c    pt. 2 is unprimed, so this is technically the az
	c
	"""
        rhs1 = (a - dd) * (a - dd) + (b - ee) * (b - ee) + c * c - 2.
        rhs2 = (a - gg) * (a - gg) + (b - hh) * (b - hh) + (c - kk) * (c - kk) - 2.
        daz = np.arctan2(rhs1, rhs2)

        daz_idx = np.atleast_1d(daz < 0.0).nonzero()[0]
        if len(daz_idx) != 0:
            if np.isscalar(daz):
                daz += 2 * np.pi
            else:
                daz[daz_idx] += 2 * np.pi

        self.az = daz / rad
        """
	c
	c   Make sure 0.0 is always 0.0, not 360.
	c
	"""
        # idx = np.where(np.abs(self.baz - 360.) < .00001)[0]
        idx = np.atleast_1d(np.abs(self.baz - 360.) < .00001).nonzero()[0]
        if len(idx) != 0:
            if np.isscalar(self.baz):
                self.baz = 0.0
            else:
                self.baz[idx] = 0.0
        # idx = np.where(np.abs(self.baz) < .00001)[0]
        idx = np.atleast_1d(np.abs(self.baz) < .00001).nonzero()[0]
        if len(idx) != 0:
            if np.isscalar(self.baz):
                self.baz = 0.0
            else:
                self.baz[idx] = 0.0

        # idx = np.where(np.abs(self.az - 360.) < .00001)[0]
        idx = np.atleast_1d(np.abs(self.az - 360.) < .00001).nonzero()[0]
        if len(idx) != 0:
            if isinstance(self.az, float):
                self.az = 0.0
            else:
                self.az[idx] = 0.0
        # idx = np.where(np.abs(self.az) < .00001)[0]
        idx = np.atleast_1d(np.abs(self.az) < .00001).nonzero()[0]
        if len(idx) != 0:
            if isinstance(self.az, float):
                self.az = 0.0
            else:
                self.az[idx] = 0.0
        
        # la_idx = np.where(lat1 == lat2)[0]
        # lo_idx = np.where(lon1 == lon2)[0]
        la_idx = np.atleast_1d(lat1 == lat2).nonzero()[0]
        lo_idx = np.atleast_1d(lon1 == lon2).nonzero()[0]
        idx = np.intersect1d(la_idx, lo_idx)
        if len(idx) != 0:
            if isinstance(self.delta, float):
                self.delta = 0.
            else:
                self.delta[idx] = 0.
            if isinstance(self.az, float):
                self.az = 0.
            else:
                self.az[idx] = 0.
            if isinstance(self.baz, float):
                self.baz = 0.
            else:
                self.baz[idx] = 0.

    def getDelta(self):
        return self.delta

    def getAz(self):
        return self.az

    def getBaz(self):
        return self.baz

    def degreesToKilometers(self):
        return self.delta * 111.19
