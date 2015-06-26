from datetime import datetime, timedelta
from math import (floor, pi, radians, degrees, sin, cos,
                  acos, sqrt, copysign, atan, asin, atan2)
import requests
import sys
import time
import pickle

def ymd2jd(year, month, day):
    # need to find the source for this
    if month == 1 or month == 2:
        yprime = year - 1
        mprime = month + 12
    else:
        yprime = year
        mprime = month

    if year > 1582 or (year == 1582 and (month >= 10 and day >= 15)):
        A = int(yprime / 100)
        B = 2 - A + int(A/4.0)
    else:
        B = 0

    if yprime < 0:
        C = int((365.25 * yprime) - 0.75)
    else:
        C = int(365.25 * yprime)

    D = int(30.6001 * (mprime + 1))

    return B + C + D + day + 1720994.5


def utcDatetime2gmst(datetimeObj):
    # need to find the source for this
    # Returns GMST in degrees
    jd = ymd2jd(datetimeObj.year, datetimeObj.month, datetimeObj.day)

    S = jd - 2451545.0
    T = S / 36525.0
    T0 = 6.697374558 + (2400.051336 * T) + (0.000025862 * T**2)
    T0 = T0 % 24

    time = datetimeObj.time()
    hour = time.hour
    minute = time.minute/60.0
    second = time.second / 3600.0
    microsecond = time.microsecond/3600000000

    UT = (hour + minute + second + microsecond) * 1.002737909
    T0 += UT
    GST = (T0 % 24) * 15

    return GST


def downloadTLE():
    # Download TLE data from all the URL's in sitelist.txt
    # And return parsed TLE's as dict TLE

    TLE = {}
    # Load all TLE sites from file.
    try:
        siteList = open('sitelist.txt', 'r')
    except IOError:
        print 'List of TLE sites could not be opened'
        sys.exit(1)

    URLlist = siteList.read().splitlines()
    with open('TLE.pickle', 'wb') as handle:
        for URL in URLlist:
            print('Downloading from ' + URL)
            response = requests.get(URL)
            data = response.content
            TLE.update(interpretTLE(data))

        # Now add a timestamp
        TLE.update({'TIMESTAMP': datetime.utcnow().strftime("%Y %m %d %H:%M")})
        pickle.dump(TLE, handle)


def loadTLE():
    # Check if we have a TLE file that is recent, or download
    # a new one.

    try:
        handle = open('TLE.pickle', 'rb')
    except IOError:
        # TLE file missing, create it
        downloadTLE()

    handle = open('TLE.pickle', 'rb')
    # Depickle
    TLE = pickle.load(handle)

    if (datetime.utcnow() - datetime.strptime(TLE['TIMESTAMP'], '%Y %m %d %H:%M')) > timedelta(days=3):
        # Pickle is too old
        handle.close()
        downloadTLE()

    with open('TLE.pickle', 'rb') as handle:
        return pickle.load(handle)


def interpretTLE(data):
    # Take a raw set of TLE's and store them in a dict
    # Split the data into lines

    if type(data) is not list:
        data = data.split('\n')
    TLE = {}

    for i in range(0, len(data)-2, 3):
        sat_name = data[i].rstrip()
        line1 = data[i+1].split()
        line2 = data[i+2].split()

        TLE[sat_name] = {'CATALOGNUM': line1[1],
                         'EPOCHTIME': float(line1[3]),
                         'DECAY': line1[4],# FIX used to be float but new data seems to have - in it
                         'INCLINATION': float(line2[2]),
                         'RAAN': float(line2[3]),
                         'ECCENTRICITY': float(line2[4])/10000000.0,
                         'ARGPERIGEE': float(line2[5]),
                         'MNANOM': float(line2[6]),
                         'MNMOTION': float(line2[7][:10]),
                         'ORBITNUM': line2[7][11:16]}

    return TLE


def pol2cart(r, geocentricRA, geocentricDec):
    # Takes inputs in radians, returns in km
    # Will give location of satellite in X Y Z geocentric

    Xg = r * cos(geocentricRA)*cos(geocentricDec)
    Yg = r * sin(geocentricRA)*cos(geocentricDec)
    Zg = r * sin(geocentricDec)

    return [Xg, Yg, Zg]


def LLA2cart(obsvLLA):
    obsvlat, obsvlong, height = obsvLLA
    # Lat long need to be in degrees, altitude in m
    # Converts lat/long into X Y Z coords in km

    a = 6378.137
    b = 6356.75231424518
    ecc = sqrt((a**2 - b**2)/a**2)
    N = a/sqrt(1 - ecc**2 * sin(obsvlat)**2)
    X = (N + height) * cos(obsvlat) * cos(obsvlong)
    Y = (N + height) * cos(obsvlat) * sin(obsvlong)
    Z = (b**2/a**2 * N + height)*sin(obsvlat)

    return [X, Y, Z]


def cart2RADec(satPos, obsvPos):
    # Returns in radians
    xg, yg, zg = satPos
    ag, bg, cg = obsvPos

    xs = xg - ag
    ys = yg - bg
    zs = zg - cg

    if ys > 0 and xs > 0:
        alpha = atan2(ys, xs)
    elif xs < 0:
        alpha = pi + atan2(ys, xs)
    elif ys < 0 and xs > 0:
        alpha = 2*pi - atan2(ys, xs)

    r = sqrt(xs**2 + ys**2 + zs**2)
    delta = asin(zs/r)

    return [alpha, delta, r]


def LLA2AzEl(satLLACoords, obsvLLA, satPos, obsvPos):

    obsvLat, obsvLong, height = obsvLLA
    satLat, satLong = satLLACoords
    satLat = radians(satLat)
    satLong = radians(satLong)

    dlong = satLong - obsvLong
    arg1 = sin(dlong)*cos(satLat)
    arg2 = cos(obsvLat)*sin(satLat) - sin(obsvLat)*cos(satLat)*cos(dlong)
    azimuth = atan2(arg1, arg2)
    azimuth = (azimuth + 2 * pi) % (2*pi)

    xd, yd, zd = satPos
    x, y, z = obsvPos

    dx = xd - x
    dy = yd - y
    dz = zd - z

    magnitude = sqrt((x**2+y**2+z**2)*(dx**2+dy**2+dz**2))
    elevation = pi/2.0 - acos((x*dx + y*dy + z*dz) / magnitude)

    return [azimuth, elevation]


def ECEF2LLA(pos):
    X, Y, Z = pos

    r = sqrt(X**2 + Y**2 + Z**2)
    p = sqrt(X**2 + Y**2)

    # First calculate geocentric lat
    lati = atan2(p, Z)
    long = atan2(Y, X)

    error = 1111
    # Earth radius at equator and poles, metres
    a = 6378137.0
    b = 6356752.31424518
    e = sqrt((a**2 - b**2)/(a**2))
    i = 0

    while(error > 0.0000001):
        if(i > 500):
            lati = atan(Z/sqrt(X**2 + Y**2))

            while(long > pi):
                long -= 2*pi
            while(long < -pi):
                long += 2*pi
            return [degrees(lati), degrees(long)]

        Rn = a/((1-(e**2)*(sin(lati)**2)))
        h = p/cos(lati) - Rn*lati
        latnext = atan(Z/p * ((1-(e**2) * Rn/(Rn + h)))**-1)
        error = abs(lati - latnext)
        lati = latnext
        i += 1

    # Bring our longitude back within 0:2pi
    while(long > pi):
        long -= 2*pi
    while(long < -pi):
        long += 2*pi

    return [degrees(lati), degrees(long)]


def getAzEl(TLE, satName, obsvLLA):

    sat = TLE[satName]
    # FIXASA obsvCoords = LLA2cart(obsvLLA)
    satLLACoords, cartCoords, obsvCoords = TLE2coords(sat, obsvLLA)
    azimuth, elevation = LLA2AzEl(satLLACoords, obsvLLA, cartCoords, obsvCoords)

    return [degrees(azimuth), degrees(elevation)]


class Satellite(object):

    def __init__(self, name, satelliteTLE):
        self.name = name
        self.catalogNumber = satelliteTLE['CATALOGNUM']
        self.epoch = satelliteTLE['EPOCHTIME']
        self.decay = satelliteTLE['DECAY']
        self.inclination = satelliteTLE['INCLINATION']
        self.RAAN = satelliteTLE['RAAN']
        self.eccentricity = satelliteTLE['ECCENTRICITY']
        self.argPerigee = satelliteTLE['ARGPERIGEE']
        self.meanAnomaly = satelliteTLE['MNANOM']
        self.meanMotion = satelliteTLE['MNMOTION']
        self.orbitNumber = satelliteTLE['ORBITNUM']
        self.SemiMajorAxis = self.calcSemiMajorAxis()
        self.perigee = self.calcPerigee()

    def runCalcs(self, manualTime=0):
        self.time = self.epochDiff() + manualTime/86400.0
        self.meanAnomaly_t = self.calcMeanAnomaly(self.time)
        self.eccAnom = self.calcEA()
        self.trueAnomaly = self.calcTrueAnomaly()
        self.geocentricDistance = self.calcGeocentricDistance()
        self.pRAAN, self.pAP = self.calcPrecession(self.time)
        self.argLatitude = self.calcArgLat()
        self.RADiff = self.calcRADiff()
        self.geocentricRA, self.geocentricDec = self.calcGeocentricRADec(manualTime)

    def calcSemiMajorAxis(self):
        # Returns distance from earth center in km
        # mu is in km^3/day^2, gravitational coefficient
        mu = 2.97554e15
        return (mu/(2*pi*self.meanMotion)**2)**(1.0/3)

    def calcPerigee(self):
        # Returns distance from earth center in km
        return self.SemiMajorAxis*(1-self.eccentricity)

    def calcMeanAnomaly(self, timeDelta):
        # Returns a mean anomaly between 0:2*pi
        M_0 = self.meanAnomaly
        n = self.meanMotion
        t = timeDelta
        M_t = M_0 + 360*(n*t - int(n*t) - int((M_0 + 360*(n*t - int(n*t)))/360))

        return radians(M_t)

    def calcEA(self):
        # I think this returns it in radians, seems consistent at least
        E0 = self.meanAnomaly
        e = self.eccentricity
        M = self.meanAnomaly_t
        # Our error term
        E_error = 100000000
        i = 0
        while(E_error > 0.000000001):
            if i > 50:
                # If we aren't converging go for a best guess
                Edash = e**2 * sin(M)*cos(M)
                Edashdash = 1.0/2 * e**3 * sin(M) * (3*cos(M)**2 - 1)
                EccAnom = M + e*sin(M) + Edash + Edashdash
                return EccAnom
            EccAnom = E0 - (E0 - e*sin(E0) - M)/(1-e*cos(E0))
            E_error = abs(EccAnom-E0)
            E0 = EccAnom
            i += 1

        return EccAnom

    def calcTrueAnomaly(self):
        # This returns TrueAnomaly in radians
        e = self.eccentricity
        x = sqrt(1-e)*cos(self.eccAnom/2)
        y = sqrt(1+e)*sin(self.eccAnom/2)
        TrueAnomaly = 2*atan2(y, x)

        return TrueAnomaly

    def calcGeocentricDistance(self):
        # Inputs: satellite data, perigee(km), true anomaly(rad)
        e = self.eccentricity
        geocentricDistance = (self.perigee*(1+e))/(1+e*cos(self.trueAnomaly))

        return geocentricDistance

    def calcArgLat(self):
        # Inputs should be in radians, returns in radians
        lat = self.pAP + self.trueAnomaly - 2*pi*(int((self.pAP + self.trueAnomaly)/(2*pi)))

        return lat

    def calcPrecession(self, timeDelta):
        # Returns in radians
        # i = radians
        # n,e,RAAN,AP,SemiMajorAxis are in their standard units from sat

        n = self.meanMotion
        e = self.eccentricity
        i = radians(self.inclination)
        RAAN = self.RAAN
        AP = self.argPerigee

        # Radius of earth in km
        Re = 6378.135
        aE = Re
        a1 = self.SemiMajorAxis/Re
        # Second gravity harmonic, J2 Propagation
        J2 = 0.0010826267

        k2 = 1/2.0 * J2 * aE**2
        d1 = 3/2.0 * k2/self.SemiMajorAxis**2 * (3*cos(i)**2 - 1)/((1-e**2)**(1.0/3))

        a0 = Re * a1 *(1 - d1/3.0 - d1**2 - 134/81.0 * d1**3)
        p0 =  Re * a0 * (1-e**2)

        # Precession of the RAAN
        # del_pRAAN will be in decimal degrees
        del_pRAAN = 360*(-3*J2*Re**2*n*timeDelta*cos(i)/(2*p0**2))
        pRAAN = RAAN + del_pRAAN

        # Precession of arg perigee
        # del_pAP will be in decimal degrees
        del_pAP = 360*(3*J2*Re**2*n*timeDelta*(5*cos(i)**2 - 1)/(4*p0**2))
        pAP = AP + del_pAP

        return [radians(pRAAN), radians(pAP)]

    def calcRADiff(self):
        # ArgLat should be in radians, Returns in radians
        i = radians(self.inclination)

        if 0 < i < pi/2 and 0 < self.argLatitude < pi:
            return acos(cos(self.argLatitude)/sqrt(1-sin(i)**2 * sin(self.argLatitude)**2))
        elif pi/2 < i < pi and pi < self.argLatitude < 2*pi:
            return acos(cos(self.argLatitude)/sqrt(1-sin(i)**2 * sin(self.argLatitude)**2))
        return 2*pi - acos(cos(self.argLatitude)/sqrt(1-sin(i)**2 * sin(self.argLatitude)**2))

    def calcGeocentricRADec(self, manual_t=0):
        # Partially working, now accounts for rotation of earth
        # Inputs should be in radians, returns in radians
        GMST = utcDatetime2gmst(datetime.utcnow())*86400.0/86164.0
        earthRotation = radians(GMST + manual_t*360.0)
        geocentricRA = -earthRotation + self.RADiff + self.pRAAN - 2*pi*(int((self.RADiff + self.pRAAN)/(2*pi)))

        geocentricDec = (copysign(1, sin(self.argLatitude))) * acos(cos(self.argLatitude)/cos(self.RADiff))

        return [geocentricRA, geocentricDec]

    def epochDiff(self):
        # Calculate difference between current time and the
        # time the TLE was set (epoch) in solar days

        # Need to rebuild this FIX
        epoch = str(self.epoch)
        epochYear = int('20'+epoch[:2])
        epochDays = float(epoch[2:])

        date_days = int(epochDays)
        date_hours = floor((epochDays - date_days)*24)
        date_min = (epochDays - date_days - date_hours/24)*1440
        date_sec = (epochDays - date_days - date_hours/24 - date_min/1440)*86400

        epochYear = str(epochYear)
        days = str(date_days)
        hours = str(int(date_hours))
        min = str(int(date_min))
        sec = str(int(date_sec))
        date_epoch = epochYear + ' ' + days + ' ' + hours + ' ' + min + ' ' + sec

        epoch_time = datetime.strptime(date_epoch, "%Y %j %H %M %S")

        curr_time = datetime.utcnow()
        # Time delta between now and epoch is therefore, in seconds
        timedelta = (curr_time - epoch_time)

        return timedelta.total_seconds()/86400

    def LLAcoordinates(self, manualTime):
        # Takes a TLE and calculates XYZ, LLA of satellite
        self.runCalcs(manualTime) # Ensure everything is up to date

        cartesianCoords = pol2cart(self.geocentricDistance, self.geocentricRA, self.geocentricDec)
        LLACoords = ECEF2LLA(cartesianCoords)

        return [LLACoords, cartesianCoords]
