import datetime
from math import pi,degrees,sqrt,sin,cos,acos,atan,atan2,tan,radians
import threading
import pytz
import requests
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
    G = 0.00000000006673
    M = 5972000000000000000000000
    P = orbitalPeriod(sat)
    a = ((P**2) * M * G/(4*pi**2 ))**(1.0/3)/1000

    return a
    
    
def calcMeanEcc(sat):
    #I think this returns it in radians, seems consistent at least
    E0 = calcMeanAnomaly(sat)
    e = sat['ECCENTRICITY']
    
    M = sat['MNANOM']
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
    a = calcSemiMajor(sat)
    e = sat['ECCENTRICITY']
    E = calcMeanEcc(sat)
    # Radius calculations
    r = a*(1-e*cos(E))-6378
    return r
    
    
def calcVel(sat):
    # Returns velocity of satellite in km/s
    #This assumes a circular orbit
    alt = calcAltitude(sat)
    a = calcSemiMajor(sat)
    mu = 398600.4
    v = sqrt(mu*(2/(alt+6378) - 1/a))

    return v


def calcTrueAnomaly(sat):
    #DONT THINK THIS IS WORKING
    e = sat['ECCENTRICITY']
    E = calcMeanEcc(sat)
    
    if E < 180:
        v = 2*atan( ((1+e)/(1-e))**(1/2) * tan(E/2))
    return degrees(v)
    
    
def loadTLE():
    URL = 'http://www.amsat.org/amsat/ftp/keps/current/nasabare.txt'
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
           'ECCENTRICITY':float(line2[4])/10000000,
           'ARGPERIGEE':float(line2[5]),
           'MNANOM':float(line2[6]),
           'MNMOTION':float(line2[7][:10]),
           'ORBITNUM':line2[7][11:16]}
           
    return TLE
        
        
def calcAngMom(sat):
    mu = 398600.4
    
    a = calcSemiMajor(sat)
    e = sat['ECCENTRICITY']
    
    h = (mu*a*(1-e**2))**(1.0/2)
    return h
    
    
#POSITION VECTORS
def posVec(sat):
    r = calcAltitude(sat)
    ohm = sat['RAAN']
    omg = sat['ARGPERIGEE']
    i = sat['INCLINATION']
    v = calcTrueAnomaly(sat)
    
    X = r*(cos(ohm)*cos(omg + v) - sin(ohm)*sin(omg+v)*cos(i))
    Y = r*(sin(ohm)*cos(omg + v) + cos(ohm)*sin(omg+v)*cos(i))
    Z = r*(sin(i)*sin(omg+v))
    
    return [X, Y, Z]


def velVec(sat):
    posVector = posVec(sat)
    X = posVector[0]
    Y = posVector[1]
    Z = posVector[2]
    
    r = calcAltitude(sat)
    ohm = radians(sat['RAAN'])
    omg = radians(sat['ARGPERIGEE'])
    e = sat['ECCENTRICITY']
    i = radians(sat['INCLINATION'])
    v = radians(calcTrueAnomaly(sat))
    h = calcAngMom(sat)
    p = r*(1+e*cos(v))
    
    Xdash = X*h*e*sin(v)/(r*p) - h/r*(cos(ohm)*sin(omg + v) + sin(ohm)*cos(omg+v)*cos(i))
    Ydash = Y*h*e*sin(v)/(r*p) - h/r*(sin(ohm)*sin(omg + v) + cos(ohm)*cos(omg+v)*cos(i))
    Zdash = Z*h*e*sin(v)/(r*p) - h*sin(i)*cos(omg + v)/r
    
    return [Xdash, Ydash, Zdash]
#print sat
#print epochDiff(sat)
#print epochDiff(sat)
#print orbitalPeriod(sat)
#print secondsIntoOrbit(sat)
#print calcMeanAnomaly(sat)
#print calcSemiMajor(sat)
#print calcMeanEcc(sat)
#print calcAltitude(sat)
#print calcVel(sat)
TLE = loadTLE()
sat = TLE['SO-50']
#print orbitalPeriod(sat)
#print calcSemiMajor(sat)
#print calcMeanEcc(sat)
#print calcTrueAnomaly(sat)
#print calcAngMom(sat)
#print posVec(sat)
#print sat['ECCENTRICITY']
#print calcTrueAnomaly(sat)
print ('posVec',posVec(sat))
print('OrbPeriod',orbitalPeriod(sat))
print('alt', calcAltitude(sat))
print('vel', calcVel(sat))
print('calcSemiMajor',calcSemiMajor(sat))
print('secondsIntoOrbit',secondsIntoOrbit(sat))
print('calcMeanAnomaly',calcMeanAnomaly(sat))
print('calcMeanEcc',calcMeanEcc(sat))
print('calcTrueAnomaly',calcTrueAnomaly(sat))
print sat

    
