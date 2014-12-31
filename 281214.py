import datetime
from math import pi,degrees,sqrt,sin,cos,acos,atan,atan2,tan,radians,asin
import pytz
import requests
import numpy

from collections import defaultdict
from time import sleep
def epochDiff(sat):
    epoch = str(sat['EPOCHTIME'])
    epochYear = int('20'+epoch[:2])
    epochDays = float(epoch[2:])

    curr_seconds = float(datetime.datetime.now(pytz.utc).strftime('%S'))
    curr_minutes = float(datetime.datetime.now(pytz.utc).strftime('%M'))
    curr_hours = float(datetime.datetime.now(pytz.utc).strftime('%H'))
    curr_days = float(datetime.datetime.now(pytz.utc).strftime('%j'))

    curr_time =  curr_days + curr_hours/24 + curr_minutes/1440 + curr_seconds/86400 

    #Time delta between now and epoch is therefore, in seconds
    timedelta = (curr_time - epochDays) * 86400
    return timedelta

def orbitalPeriod(sat):
    # Returns the orbital period in seconds
    period = 86400 / sat['MNMOTION']
    return period


def secondsIntoOrbit(sat):
    seconds = epochDiff(sat) % orbitalPeriod(sat)
    return seconds


def calcMeanAnomaly(sat):
    M=sat['MNANOM']+360*secondsIntoOrbit(sat)/orbitalPeriod(sat)
    while(M>360):
        M -= 360    
    return M
    

def calcSemiMajor(sat):
    # Returns the semi major axis in km
    G = 6.67384E-11
    M = 5.972E24
    P = orbitalPeriod(sat)
    a = ((P**2) * M * G/(4*pi**2 ))**(1.0/3)/1000

    return a
    
    
def calcEccAnom(sat):
    #I think this returns it in radians, seems consistent at least
    E0 = calcMeanAnomaly(sat)
    e = sat['ECCENTRICITY']

    M = calcMeanAnomaly(sat)
    #Initial guess
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
    #Both give the same result pretty much
    #print a*(1-e*cos(E))-6378
    r = p/(1+e*cos(v))-6378
    return r
    
    
def calcVel(sat):
    # Returns velocity of satellite in km/s
    #This assumes a circular orbit
    alt = calcAltitude(sat)
    a = calcSemiMajor(sat)
    mu = 398600.44189
    v = sqrt(mu*(2/(alt+6378) - 1/a))

    return v


def calcTrueAnomaly(sat):
    #DONT THINK THIS IS WORKING
    e = sat['ECCENTRICITY']
    E = radians(calcEccAnom(sat))
    
    v = 2*atan( ((1+e)/(1-e))**(1/2) * tan(E/2))
    if v>pi:
        v = v - 2*pi
    return degrees(v)
    
    
def loadTLE():
    URL = 'http://www.amsat.org/amsat/ftp/keps/current/nasabare.txt'
    # URL = 'http://www.celestrak.com/NORAD/elements/stations.txt'
    response = requests.get(URL)
    data = response.content
    TLE = defaultdict()

    # split the data into lines
    data = data.split('\n')

    for i in range(0,len(data)-2,3):
        line1 = data[i+1].split()
        line2 = data[i+2].split()

        TLE[data[i]] = {'CATALOGNUM':line1[1],
           'EPOCHTIME':float(line1[3]),
           'DECAY':float(line1[4]),
           'INCLINATION':float(line2[2]),
           'RAAN':float(line2[3]),
           'ECCENTRICITY':float(line2[4])/10000000.0,
           'ARGPERIGEE':float(line2[5]),
           'MNANOM':float(line2[6]),
           'MNMOTION':float(line2[7][:10]),
           'ORBITNUM':line2[7][11:16]}
           
           
    TLE['TEST'] = {'CATALOGNUM':'88888U',
           'EPOCHTIME':14364.07403355,
           'DECAY':.00013787,
           'INCLINATION':51.6453,
           'RAAN':208.0014,
           'ECCENTRICITY':0.0007496,
           'ARGPERIGEE':187.5035,
           'MNANOM':321.7881,
           'MNMOTION':15.52846000,
           'ORBITNUM':999}
    return TLE
        
        
def calcAngMom(sat):
    mu = 398600.44189
    a = calcSemiMajor(sat)
    e = sat['ECCENTRICITY']
    
    h = (mu*a*(1-e**2))**(1.0/2)
    return h
    
    
#POSITION VECTORS
def posVec(sat):
    #DEPRECATED, GOING TO BE REMOVED
    r = calcAltitude(sat)
    ohm = radians(sat['RAAN'])
    omg = radians(sat['ARGPERIGEE'])
    i = radians(sat['INCLINATION'])
    v = radians(calcTrueAnomaly(sat))
    
    X = r*(cos(ohm)*cos(omg + v) - sin(ohm)*sin(omg+v)*cos(i))
    Y = r*(sin(ohm)*cos(omg + v) + cos(ohm)*sin(omg+v)*cos(i))
    Z = r*(sin(i)*sin(omg+v))
    
    return [X, Y, Z]


def J2000():
    #DEPRECATED, GOING TO BE REMOVED
    d0 = 2451545.0
    d1 = datetime.datetime.utcfromtimestamp(0)
    print d1
    #delta = d1 - d0
    #return delta
    

def velVec(sat):
    #DEPRECATED, GOING TO BE REMOVED
    posVector = posVec(sat)
    X = posVector[0]
    Y = posVector[1]
    Z = posVector[2]
    mu = 398600.44189
    
    r = calcAltitude(sat)
    ohm = radians(sat['RAAN'])
    omg = radians(sat['ARGPERIGEE'])
    e = sat['ECCENTRICITY']
    i = radians(sat['INCLINATION'])
    v = radians(calcTrueAnomaly(sat))
    h = calcAngMom(sat)
    p = h**2/mu
    
    Xdash = X*h*e*sin(v)/(r*p) - h/r*(cos(ohm)*sin(omg + v) + sin(ohm)*cos(omg+v)*cos(i))
    Ydash = Y*h*e*sin(v)/(r*p) - h/r*(sin(ohm)*sin(omg + v) + cos(ohm)*cos(omg+v)*cos(i))
    Zdash = Z*h*e*sin(v)/(r*p) - h*sin(i)*cos(omg + v)/r
    
    return [Xdash, Ydash, Zdash]

    
def ECItoECF(pos,vel):
    #DEPRECATED, GOING TO BE REMOVED
    X = pos[0]
    Y = pos[1]
    Z = pos[2]
    Xdash = vel[0]
    Ydash = vel[1]
    Zdash = vel[2]
    
    
    # is the Greenwich hour angle of the Earth's prime meridian,
    D0 = 2457022.125
    GST = 18.697374558 + 24.06570982441908 * D0
    while(GST>360):
        GST -=360
    theta = radians(GST)
    #Rotation rate of earth
    w = 7.2921158553E-5 

    Xecf = cos(theta)*(Xdash + w*Y) + sin(theta)*(Ydash - w*X)
    Yecf = -1*sin(theta)*(Xdash + w*Y) + cos(theta)*(Ydash - w*X)
    Zecf = Z
    
    ECF = [Xecf, Yecf, Zecf]
    return ECF
    

def ECEF2LLA(pos):
    #Non working lat, long works
    X = pos[0]
    Y = pos[1]
    Z = pos[2]
    a = 6378137
    e = 8.1819190842622e-2
    b = (a**2 * (1-e**2))**(1/2)
    
    r = calcAltitude(sat)
    ep = ((a**2 - b**2)/(b**2))**(1/2)
    
    p = (X**2 + Y**2)**(1/2)
    th = atan2(a*Z,b*p)

    
    
    #2.2EceFtoLLA microem.ru, iterative method
    long = atan2(Y,X)
    
    H0 = 0
    lati = atan(Z/(p*(1-e**2))) #Consider atan2
    
    error = 100
    while(error>0.00001):
        Ni = a/(1 - (e**2)*(sin(lati)**2))**(1/2)
        hip1 = p/cos(lati) - Ni
        latip1 = atan(Z/(p*(1 - (e**2)*Ni/(Ni + hip1))))
        error = abs(lati - latip1)
        lati = latip1
        print degrees(lati)
    
    lat = lati
    
    
    
    a = 6378137
    b = 6356752.31424518
    e = ((a**2 - b**2)/a**2)**(1/2)
    edash = ((a**2 - b**2)/b**2)**(1/2)
    theta = atan2(Z*a,(p*b))
    #Try closed form instead
    #lat = atan2((Z + (edash**2)*b*(sin(theta))**3),(p - (e**2)*a*(cos(theta)**3)))
    
    out = [degrees(lat), degrees(long)]
    return out
    
    
def LLA2ECEF(lat,lon,alt):
    #CONFIRMED WORKING
    # see http://www.mathworks.de/help/toolbox/aeroblks/llatoecefposition.html
    lat = radians(lat)
    lon = radians(lon)
    radius = 6378137.0 #radius of earth in major
    #flattening factor
    f = 1.0/298.257223563
    
    cosLat = cos(lat)
    sinLat = sin(lat)
    FF = (1.0 - f)**2
    C = 1/sqrt(cosLat**2 + FF * sinLat**2)
    S = C * FF
    
    X = (radius * C + alt)* cosLat * cos(lon)
    Y = (radius * C + alt)* cosLat * sin(lon)
    Z = (radius * S + alt)* sinLat
    
    return [X/1000, Y/1000, Z/1000]
    
def calcAzi(satpos,recvpos):
    X = satpos[0] - recvpos[0]
    Y = satpos[1] - recvpos[1]
    angle = degrees(atan(Y/X))
    return angle
    
def TLEtoVec(sat):
    mu = 398600.44189
    a = calcSemiMajor(sat)
    e = sat['ECCENTRICITY']
    v = radians(calcTrueAnomaly(sat))
    AP = radians(sat['ARGPERIGEE'])
    inc = radians(sat['INCLINATION'])
    RAAN = radians(sat['RAAN'])
    
    
    slr = a*(1-e**2)
    rm = slr / (1 + e * cos(v))
    
    arglat = AP + v
    
    sarglat = sin(arglat)
    carglat = cos(arglat)
    
    c4 = (mu/slr)**(1/2)
    c5 = e * cos(AP) + carglat
    c6 = e * sin(AP) + sarglat
    
    sinc = sin(inc)
    cinc = cos(inc)
    
    sRAAN = sin(RAAN)
    cRAAN = cos(RAAN)
    
    #Position vector
    X = rm * (cRAAN * carglat - sRAAN * cinc * sarglat)
    Y = rm * (sRAAN * carglat + cRAAN * cinc * sarglat)
    Z = rm * sinc * sarglat
    
    # Velocity vector
    Xdash = -c4 * (cRAAN * c6 + sRAAN * cinc * c5);
    Ydash = -c4 * (sRAAN * c6 - cRAAN * cinc * c5);
    Zdash = c4 * c5 * sinc;
    
    pos = [X, Y, Z]
    vel = [Xdash, Ydash, Zdash]
    return [pos, vel]



TLE = loadTLE()
TLE['TEST']['customtime']=1
sat = TLE['SO-50']



if(2==1):
    print('posVec',posVec(sat))
    print('OrbPeriod',orbitalPeriod(sat))
    print('alt', calcAltitude(sat))
    print('vel', calcVel(sat))
    print('calcSemiMajor',calcSemiMajor(sat))
    print('secondsIntoOrbit',secondsIntoOrbit(sat))
    print('calcMeanAnomaly',calcMeanAnomaly(sat))
    print('calcMeanEcc',calcEccAnom(sat))
    print('calcTrueAnomaly',radians(calcTrueAnomaly(sat)))
    print('EpochDiff',epochDiff(sat))
    print('RAAN',radians(sat['RAAN']))
    print('Eccentricity',sat['ECCENTRICITY'])
    print('Inclination',radians(sat['INCLINATION']))
    print('Argperigee',radians(sat['ARGPERIGEE']))


testfile = open('output.txt','w')
i = 0;
while(i+10000000<orbitalPeriod(sat)/60):

    v = calcTrueAnomaly(sat)
    r = calcAltitude(sat)
    x = r * cos(v)
    y = r * sin(v)
    pas = posVec(sat)
    data = "{} {} {};\n".format(pas[0],pas[1],pas[2])
    testfile.write(data)
    print(i)
    i += 1

    
satpos = TLEtoVec(sat)[0]
recvpos = LLA2ECEF(30.331446,145.406932,112)
#print calcAzi(satpos, recvpos)
print recvpos
print ECEF2LLA(recvpos)