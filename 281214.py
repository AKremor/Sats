import datetime
from math import pi, degrees, sqrt, sin, cos, acos, atan, atan2, tan, radians
import pytz
import requests
from time import sleep
from numpy import (sin,cos,tan,sqrt,radians,arctan2,hypot,degrees,mod,
                   atleast_2d,atleast_1d,empty_like,asarray,array)

from random import random
def epochDiff(sat):
    epoch = str(sat['EPOCHTIME'])
    epochYear = int('20'+epoch[:2])
    epochDays = float(epoch[2:])

    curr_seconds = float(datetime.datetime.now(pytz.utc).strftime('%S'))
    curr_minutes = float(datetime.datetime.now(pytz.utc).strftime('%M'))
    curr_hours = float(datetime.datetime.now(pytz.utc).strftime('%H'))
    curr_days = float(datetime.datetime.now(pytz.utc).strftime('%j'))

    curr_time = curr_days + curr_hours/24 + curr_minutes/1440 + curr_seconds/86400

    # Time delta between now and epoch is therefore, in seconds
    timedelta = (curr_time - epochDays) * 86400
    return timedelta


def orbitalPeriod(sat):
    # Returns the orbital period in seconds
    period = 86400 / sat['MNMOTION']
    return period


def secondsIntoOrbit(sat):
    # Returns how many seconds into orbit sat is
    seconds = epochDiff(sat) % orbitalPeriod(sat)
    return seconds


def calcMeanAnomaly(sat):
    M = sat['MNANOM'] + 360 * secondsIntoOrbit(sat) / orbitalPeriod(sat)
    while(M > 360):
        M -= 360
    return M


def calcSemiMajor(sat):
    # Returns the semi major axis in km
    G = 6.67384E-11
    M = 5.972E24
    P = orbitalPeriod(sat)
    a = ((P**2) * M * G/(4*pi**2))**(1.0/3)/1000

    return a


def calcEccAnom(sat):
    # I think this returns it in radians, seems consistent at least
    E0 = calcMeanAnomaly(sat)
    e = sat['ECCENTRICITY']

    M = calcMeanAnomaly(sat)
    # Initial guess
    Ed = (M + e*sin(E0) - E0)/(1-e*cos(E0))

    # Our delta E term
    Ed = 100000000
    while(Ed > 0000000000000.000000001):
        E = E0 - (E0 - e*sin(E0) - M)/(1-e*cos(E0))
        Ed = abs(E-E0)
        E0 = E

    return E


def calcAltitude(sat):
    # Returns altitude in km
    mu = 398600.44189
    e = sat['ECCENTRICITY']
    v = radians(calcTrueAnomaly(sat))
    p = calcAngMom(sat)**2/mu
    # Radius calculations
    # Both give the same result pretty much
    # print a*(1-e*cos(E))-6378
    r = p/(1+e*cos(v))-6378
    return r


def calcVel(sat):
    # Returns velocity of satellite in km/s
    # This assumes a circular orbit
    alt = calcAltitude(sat)
    a = calcSemiMajor(sat)
    mu = 398600.44189
    v = sqrt(mu*(2/(alt+6378) - 1/a))

    return v


def calcTrueAnomaly(sat):
    # DONT THINK THIS IS WORKING
    e = sat['ECCENTRICITY']
    E = radians(calcEccAnom(sat))
   
    #v = 2*atan(((1+e)/(1-e))**(1.0/2) * tan(E/2))
    v = acos((cos(E) - e)/(1 - e*cos(E)))
    
    v = 2*atan2(sqrt(1-e)*cos(E/2),sqrt(1+e)*sin(E/2))
    if v < 0:
        v = v + 2 * pi
    
    print('vnew',degrees(v))
    return degrees(v)


def loadTLE():
    URL = 'http://www.amsat.org/amsat/ftp/keps/current/nasabare.txt'
    # URL = 'http://www.celestrak.com/NORAD/elements/stations.txt'
    response = requests.get(URL)
    data = response.content
    TLE = {}

    # split the data into lines
    data = data.split('\n')

    for i in range(0, len(data)-2, 3):
        line1 = data[i+1].split()
        line2 = data[i+2].split()

        TLE[data[i]] = {'CATALOGNUM': line1[1],
                        'EPOCHTIME': float(line1[3]),
                        'DECAY': float(line1[4]),
                        'INCLINATION': float(line2[2]),
                        'RAAN': float(line2[3]),
                        'ECCENTRICITY': float(line2[4])/10000000.0,
                        'ARGPERIGEE': float(line2[5]),
                        'MNANOM': float(line2[6]),
                        'MNMOTION': float(line2[7][:10]),
                        'ORBITNUM': line2[7][11:16]}

    return TLE


def calcAngMom(sat):
    mu = 398600.44189
    a = calcSemiMajor(sat)
    e = sat['ECCENTRICITY']

    h = (mu*a*(1-e**2))**(1.0/2)
    return h

    
def ECEF2LLA(pos):
    # Confirmed working
    X = pos[0]
    Y = pos[1]
    Z = pos[2]
    
    long = atan2(Y,X)
    
    
    r = (X**2 + Y**2 + Z**2)**(1.0/2)
    p = (X**2 + Y**2)**(1.0/2)
    #First calculate geocentric lat
    lati = atan2(p,Z)
    error = 1111
    a = 6378137.0
    b = 6356752.31424518
    e = ((a**2 - b**2)/(a**2))**(1.0/2)
    i = 0
    hap = 0
    while(error > 0.0001):
        if(i>20 and hap==0):
            lati = random()
            print("enabled newton raphson modified guess")
            hap  = 1
            
        Rn = a/((1-(e**2)*(sin(lati)**2)))
        h = p/cos(lati) - Rn*lati
        latnext = atan(Z/p * ((1-(e**2) * Rn/(Rn + h)))**-1)
        error = abs(lati - latnext)
        lati = latnext
        print lati
        i += 1
    
    return [degrees(lati), degrees(long)]


def LLA2ECEF(lat, lon, alt):
    # CONFIRMED WORKING
    # see http://www.mathworks.de/help/toolbox/aeroblks/llatoecefposition.html
    lat = radians(lat)
    lon = radians(lon)
    # radius of earth in major
    radius = 6378137.0
    # flattening factor
    f = 1.0/298.257223563

    cosLat = cos(lat)
    sinLat = sin(lat)
    FF = (1.0 - f)**2
    C = 1/sqrt(cosLat**2 + FF * sinLat**2)
    S = C * FF

    X = (radius * C + alt) * cosLat * cos(lon)
    Y = (radius * C + alt) * cosLat * sin(lon)
    Z = (radius * S + alt) * sinLat

    return [X/1000, Y/1000, Z/1000]


def calcAzi(satpos, recvpos):
    X = satpos[0] - recvpos[0]
    Y = satpos[1] - recvpos[1]
    angle = degrees(atan(Y/X))
    return angle
    

def J2000():
    yeard = int(datetime.datetime.now(pytz.utc).strftime('%Y')) - 2000
    dayd = int(datetime.datetime.now(pytz.utc).strftime('%j'))
    hourd = int(datetime.datetime.now(pytz.utc).strftime('%H'))
    mind = int(datetime.datetime.now(pytz.utc).strftime('%M'))
    secd = int(datetime.datetime.now(pytz.utc).strftime('%S'))

    # Less 0.5 because epoch is at midday
    days = yeard*365.25 + dayd + hourd/24.0 + mind/1440.0 + secd/86400.0 + 0.5
    print days
    return days
    
def rottrip(ang):
    #ang = ang.squeeze()
    """ported from:
    https://github.com/dinkelk/astrodynamics/blob/master/rot3.m
    """
    return array([[cos(ang),  sin(ang), 0],
                 [-sin(ang), cos(ang), 0],
                 [0,         0,        1]])
    
def eci2ecef(eci):
    eci = asarray([eci])
    d = J2000()
    GMST = 18.697374558 + 24.065709824419*d
    print GMST
    while(GMST>24):
        GMST -= 24

    lst = [2*pi*GMST/24]
    N,trip = eci.shape
    """ported from:
    https://github.com/dinkelk/astrodynamics/blob/master/rot3.m
    """
    ecef = empty_like(eci)
    for i in range(N):
        ecef[i,:] = rottrip(lst[i]).dot(eci[i,:])
    return ecef.tolist()
    
    
def ECItoECF(pos,vel):
    # KEEP ME
    X = pos[0]
    Y = pos[1]
    Z = pos[2]
    Xdash = vel[0]
    Ydash = vel[1]
    Zdash = vel[2]
    
    
    # is the Greenwich hour angle of the Earth's prime meridian,
    d = 5478.08333
    GMST = 18.697374558 + 24.065709824419*d
    while(GMST>24):
        GMST -= 24
    theta = 2*pi*GMST/24
    #Rotation rate of earth
    w = 7.2921158553E-5 

    #Xecf = cos(theta)*(Xdash + w*Y) + sin(theta)*(Ydash - w*X)
    #Yecf = -1*sin(theta)*(Xdash + w*Y) + cos(theta)*(Ydash - w*X)
    #Zecf = Z
    i = radians(sat['INCLINATION'])
    Yecf = X * sin(theta) - Y*cos(i)*cos(theta)
    Xecf = X * cos(theta) - Y*cos(i)*sin(theta)
    Zecf = Y*sin(i)
    
    
    ECF = [Xecf, Yecf, Zecf]
    return ECF
    

def TLEtoVec(sat):
    mu = 398600.44189
    a = calcSemiMajor(sat)
    e = sat['ECCENTRICITY']
    v = radians(calcTrueAnomaly(sat))
    #v = radians(170.5994401225295)
    AP = radians(sat['ARGPERIGEE'])
    inc = radians(sat['INCLINATION'])
    RAAN = radians(sat['RAAN'])

    slr = a*(1-e**2)
    rm = slr / (1 + e * cos(v))

    arglat = AP + v

    sarglat = sin(arglat)
    carglat = cos(arglat)

    c4 = (mu/slr)**(1.0/2)
    c5 = e * cos(AP) + carglat
    c6 = e * sin(AP) + sarglat

    sinc = sin(inc)
    cinc = cos(inc)

    sRAAN = sin(RAAN)
    cRAAN = cos(RAAN)

    # Position vector
    X = rm * (cRAAN * carglat - sRAAN * cinc * sarglat)
    Y = rm * (sRAAN * carglat + cRAAN * cinc * sarglat)
    Z = rm * sinc * sarglat

    # Velocity vector
    Xdash = -c4 * (cRAAN * c6 + sRAAN * cinc * c5)
    Ydash = -c4 * (sRAAN * c6 - cRAAN * cinc * c5)
    Zdash = c4 * c5 * sinc

    pos = [X, Y, Z]
    vel = [Xdash, Ydash, Zdash]
    return [pos, vel]


TLE = loadTLE()
sat = TLE['SO-50']


if(1 == 2):
    print('posVec', TLEtoVec(sat))
    print('OrbPeriod', orbitalPeriod(sat))
    print('alt', calcAltitude(sat))
    print('vel', calcVel(sat))
    print('calcSemiMajor', calcSemiMajor(sat))
    print('secondsIntoOrbit', secondsIntoOrbit(sat))
    print('calcMeanAnomaly', calcMeanAnomaly(sat))
    print('calcMeanEcc', calcEccAnom(sat))
    print('calcTrueAnomaly', (calcTrueAnomaly(sat)))
    print('EpochDiff', epochDiff(sat))
    #print('RAAN', radians(sat['RAAN']))
    #print('Eccentricity', sat['ECCENTRICITY'])
    #print('Inclination', radians(sat['INCLINATION']))
    #print('Argperigee', radians(sat['ARGPERIGEE']))


testfile = open('output.txt', 'w')
i = 0
while(i + 10000000 < orbitalPeriod(sat)/60):

    v = calcTrueAnomaly(sat)
    r = calcAltitude(sat)
    x = r * cos(v)
    y = r * sin(v)
    pas = posVec(sat)
    data = "{} {} {};\n".format(pas[0], pas[1], pas[2])
    testfile.write(data)
    print(i)
    i += 1


#J2000()
satvec = TLEtoVec(sat)
print('satpos',satvec)
converted = eci2ecef(satvec[0])
print('ECItoECF',converted)

print('ECEF2LLA',ECEF2LLA(converted[0]))
print('calcTrueAnomaly', (calcTrueAnomaly(sat)))


# Only latitude coordsis correct,longitude conver wrong
# Lat: coords correct, my conversion wrong
# Long: coords wrong, my conversion corect