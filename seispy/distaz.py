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
    """

    def __init__(self, lat1, lon1, lat2, lon2):

        self.stalat = lat1
        self.stalon = lon1
        self.evtlat = lat2
        self.evtlon = lon2
        if (lat1 == lat2) and (lon1 == lon2):
            self.delta = 0.0
            self.az = 0.0
            self.baz = 0.0
            return

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

        scolat = math.pi / 2.0 - math.atan((1. - sph) * (1. - sph) * math.tan(lat1 * rad))
        ecolat = math.pi / 2.0 - math.atan((1. - sph) * (1. - sph) * math.tan(lat2 * rad))
        slon = lon1 * rad
        elon = lon2 * rad
        """
	c
	c  a - e are as defined by Bullen (pg. 154, Sec 10.2)
	c     These are defined for the pt. 1
	c
        """
        a = math.sin(scolat) * math.cos(slon)
        b = math.sin(scolat) * math.sin(slon)
        c = math.cos(scolat)
        d = math.sin(slon)
        e = -math.cos(slon)
        g = -c * e
        h = c * d
        k = -math.sin(scolat)
        """
	c
	c  aa - ee are the same as a - e, except for pt. 2
	c
        """
        aa = math.sin(ecolat) * math.cos(elon)
        bb = math.sin(ecolat) * math.sin(elon)
        cc = math.cos(ecolat)
        dd = math.sin(elon)
        ee = -math.cos(elon)
        gg = -cc * ee
        hh = cc * dd
        kk = -math.sin(ecolat)
        """
	c
	c  Bullen, Sec 10.2, eqn. 4
	c
        """
        delrad = math.acos(a * aa + b * bb + c * cc)
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
        dbaz = math.atan2(rhs1, rhs2)
        if (dbaz < 0.0):
            dbaz = dbaz + 2 * math.pi

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
        daz = math.atan2(rhs1, rhs2)
        if daz < 0.0:
            daz = daz + 2 * math.pi

        self.az = daz / rad
        """
	c
	c   Make sure 0.0 is always 0.0, not 360.
	c
	"""
        if (abs(self.baz - 360.) < .00001):
            self.baz = 0.0
        if (abs(self.az - 360.) < .00001):
            self.az = 0.0

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
