import math
import numpy as np


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


def km2deg(kilometers):
    return kilometers / 111.19


def deg2km(degree):
    return degree * 111.19


class distaz:
    """
    c Subroutine to calculate the Great Circle Arc distance
    c    between two sets of geographic coordinates
    c
    c Equations take from Bullen, pages 154, 155
    c
    c T. Owens, September 19, 1991
    c           Sept. 25 -- fixed az and baz calculations
    c
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

    Mijian Xu, Jan 01, 2016
    Add np.ndarray to availible input varible.
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

        rad = 2. * math.pi / 360.0
        """
        c
        c scolat and ecolat are the geocentric colatitudes
        c as defined by Richter (pg. 318)
        c
        c Earth Flattening of 1/298.257 take from Bott (pg. 3)
        c
        """
        sph = 1.0 / 298.257

        scolat = math.pi / 2.0 - np.arctan((1. - sph) * (1. - sph) * np.tan(lat1 * rad))
        ecolat = math.pi / 2.0 - np.arctan((1. - sph) * (1. - sph) * np.tan(lat2 * rad))
        slon = lon1 * rad
        elon = lon2 * rad
        """
	c
	c  a - e are as defined by Bullen (pg. 154, Sec 10.2)
	c     These are defined for the pt. 1
	c
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
        delrad = np.arccos(a * aa + b * bb + c * cc)
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

        dbaz_idx = np.where(dbaz < 0.0)[0]
        if len(dbaz_idx) != 0:
            if isinstance(dbaz, (int, float)):
                dbaz += 2 * math.pi
            else:
                dbaz[dbaz_idx] += 2 * math.pi

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

        daz_idx = np.where(daz < 0.0)[0]
        if len(daz_idx) != 0:
            if isinstance(daz, (int, float)):
                daz += 2 * math.pi
            else:
                daz[daz_idx] += 2 * math.pi

        self.az = daz / rad
        """
	c
	c   Make sure 0.0 is always 0.0, not 360.
	c
	"""
        idx = np.where(np.abs(self.baz - 360.) < .00001)[0]
        if len(idx) != 0:
            if isinstance(self.baz, float):
                self.baz = 0.0
            else:
                self.baz[idx] = 0.0
        idx = np.where(np.abs(self.baz) < .00001)[0]
        if len(idx) != 0:
            if isinstance(self.baz, float):
                self.baz = 0.0
            else:
                self.baz[idx] = 0.0

        idx = np.where(np.abs(self.az - 360.) < .00001)[0]
        if len(idx) != 0:
            if isinstance(self.az, float):
                self.az = 0.0
            else:
                self.az[idx] = 0.0
        idx = np.where(np.abs(self.az) < .00001)[0]
        if len(idx) != 0:
            if isinstance(self.az, float):
                self.az = 0.0
            else:
                self.az[idx] = 0.0
        
        la_idx = np.where(lat1 == lat2)[0]
        lo_idx = np.where(lon1 == lon2)[0]
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

# distaz = DistAz(0, 0, 1,1)
# print "%f  %f  %f" % (distaz.getDelta(), distaz.getAz(), distaz.getBaz())
if __name__ == '__main__':
    ela = 1
    elo = 1
    sla = 2
    slo = 1
    da = distaz(ela, elo, sla, slo)
    print(da.baz)
