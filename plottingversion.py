from datetime import datetime, timedelta
from math import (floor, pi, radians, degrees, sin, cos,
                  acos, sqrt, copysign, atan, asin, atan2)
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import requests
from time import sleep
import numpy as np
from selenium import webdriver


def ymd2jd(year, month, day):
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
    # Confirmed working
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

    GST = T0 % 24
    GST = GST * 15

    return GST


def loadTLE():
    # Confirmed working
    URLlist = ['http://www.amsat.org/amsat/ftp/keps/current/nasabare.txt',
               'http://www.celestrak.com/NORAD/elements/stations.txt',
               'http://www.celestrak.com/NORAD/elements/sbas.txt',
               'http://www.celestrak.com/NORAD/elements/tle-new.txt',
               'http://www.celestrak.com/NORAD/elements/military.txt',
               'http://www.celestrak.com/NORAD/elements/cubesat.txt',
               'http://www.celestrak.com/NORAD/elements/other.txt',
               'http://www.celestrak.com/NORAD/elements/education.txt',
               'http://www.celestrak.com/NORAD/elements/radar.txt',
               'http://www.celestrak.com/NORAD/elements/geodetic.txt',
               'http://www.celestrak.com/NORAD/elements/engineering.txt',
               'http://www.celestrak.com/NORAD/elements/science.txt',
               'http://www.celestrak.com/NORAD/elements/musson.txt',
               'http://www.celestrak.com/NORAD/elements/nnss.txt',
               'http://www.celestrak.com/NORAD/elements/galileo.txt',
               'http://www.celestrak.com/NORAD/elements/beidou.txt',
               'http://www.celestrak.com/NORAD/elements/gps-ops.txt',
               'http://www.celestrak.com/NORAD/elements/glo-ops.txt',
               'http://www.celestrak.com/NORAD/elements/amateur.txt',
               'http://www.celestrak.com/NORAD/elements/x-comm.txt',
               'http://www.celestrak.com/NORAD/elements/other-comm.txt',
               'http://www.celestrak.com/NORAD/elements/iridium.txt',
               'http://www.celestrak.com/NORAD/elements/orbcomm.txt',
               'http://www.celestrak.com/NORAD/elements/globalstar.txt',
               'http://www.celestrak.com/NORAD/elements/molniya.txt',
               'http://www.celestrak.com/NORAD/elements/raduga.txt',
               'http://www.celestrak.com/NORAD/elements/gorizont.txt',
               'http://www.celestrak.com/NORAD/elements/intelsat.txt',
               'http://www.celestrak.com/NORAD/elements/geo.txt',
               'http://www.celestrak.com/NORAD/elements/argos.txt',
               'http://www.celestrak.com/NORAD/elements/tdrss.txt',
               'http://www.celestrak.com/NORAD/elements/sarsat.txt',
               'http://www.celestrak.com/NORAD/elements/dmc.txt',
               'http://www.celestrak.com/NORAD/elements/resource.txt',
               'http://www.celestrak.com/NORAD/elements/noaa.txt',
               'http://www.celestrak.com/NORAD/elements/goes.txt',
               'http://www.celestrak.com/NORAD/elements/weather.txt'
               ]
    TLE = {}

    try:
        TLEfile = open('TLE.txt', 'r')
        curr_time = datetime.utcnow()
        lines = TLEfile.readlines()
        TLE_time = datetime.strptime(lines[0].rstrip(), "%Y %m %d %H:%M")
        if (curr_time - TLE_time) > timedelta(days=3):
            TLEfile.close()
            # Open instead for writing
            TLEfile = open('TLE.txt', 'w')
            TLEfile.write(datetime.utcnow().strftime("%Y %m %d %H:%M") + '\n')
            print('Current TLE older than 3 days, downloading new')
            for URL in URLlist:
                response = requests.get(URL)
                data = response.content
                TLEfile.write(data)
                TLE.update(interpretTLE(data))

            return TLE
        else:
            print('Using cached TLE: ' + str(lines[0]))
            TLE.update(interpretTLE(lines[1:]))
            return TLE
    except:
        print('Downloading TLE')
        TLEfile = open('TLE.txt', 'w')
        TLEfile.write(datetime.utcnow().strftime("%Y %m %d %H:%M" + '\n'))

        for URL in URLlist:
            response = requests.get(URL)
            data = response.content
            TLEfile.write(data)
            TLE.update(interpretTLE(data))

        return TLE


def interpretTLE(data):
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
                         'DECAY': float(line1[4]),
                         'INCLINATION': float(line2[2]),
                         'RAAN': float(line2[3]),
                         'ECCENTRICITY': float(line2[4])/10000000.0,
                         'ARGPERIGEE': float(line2[5]),
                         'MNANOM': float(line2[6]),
                         'MNMOTION': float(line2[7][:10]),
                         'ORBITNUM': line2[7][11:16]}

    return TLE


def epochDiff(sat):
    # Confirmed working
    # This function now returns in solar days
    epoch = str(sat['EPOCHTIME'])
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


def calcMA(sat, t):
    # Confirmed working
    # Returns a mean anomaly between 0:2*pi

    M0 = sat['MNANOM']
    n = sat['MNMOTION']

    Mt = M0 + 360*(n*t - int(n*t) - int((M0 + 360*(n*t - int(n*t)))/(360)))
    return radians(Mt)


def calcEA(sat, M):
    # I think this returns it in radians, seems consistent at least
    E0 = sat['MNANOM']
    e = sat['ECCENTRICITY']

    # Initial guess
    Ed = (M + e*sin(E0) - E0)/(1-e*cos(E0))

    # Our delta E term
    Ed = 100000000
    i = 0
    while(Ed > 0.000000001):
        if i > 50:
            Edash = e**2 * sin(M)*cos(M)
            Edashdash = 1.0/2 * e**3 * sin(M) * (3*cos(M)**2 - 1)
            E = M + e*sin(M) + Edash + Edashdash
            return E
        E = E0 - (E0 - e*sin(E0) - M)/(1-e*cos(E0))
        Ed = abs(E-E0)
        E0 = E
        i += 1
    return E


def calcTA(sat, M):
    # Confirmed working
    # This returns v in radians
    e = sat['ECCENTRICITY']

    E = calcEA(sat, M)
    x = sqrt(1-e)*cos(E/2)
    y = sqrt(1+e)*sin(E/2)
    v = 2*atan2(y, x)
    return v


def calcSMA(sat):
    # Confirmed working
    # Returns distance from earth center in km

    n = sat['MNMOTION']
    # mu is in km^3/day^2
    mu = 2.97554e15
    a = (mu/(2*pi*n)**2)**(1.0/3)
    return a


def calcPerigee(sat, a):
    # Confirmed working
    # Returns distance from earth center in km
    e = sat['ECCENTRICITY']

    perigee = a*(1-e)
    return perigee


def geoDist(sat, P, v):
    # Confirmed working
    e = sat['ECCENTRICITY']

    r = (P*(1+e))/(1+e*cos(v))
    return r


def precession(sat, a, t):
    # Confirmed working, returns in radians
    # i = radians
    # n,e,RAAN,AP,a are in their standard units from sat

    n = sat['MNMOTION']
    e = sat['ECCENTRICITY']
    i = radians(sat['INCLINATION'])
    RAAN = sat['RAAN']
    AP = sat['ARGPERIGEE']

    # Radius of earth in km
    Re = 6378.135
    aE = Re
    a1 = a/Re
    # Second gravity harmonic, J2 Propagation
    J2 = 0.0010826267

    k2 = 1/2.0 * J2 * aE**2
    d1 = 3/2.0 * k2/a**2 * (3*cos(i)**2 - 1)/((1-e**2)**(1.0/3))

    a0 = a1*(1 - d1/3.0 - d1**2 - 134/81.0 * d1**3)
    p0 = a0 * (1-e**2)
    a0 = Re * a0
    p0 = Re * p0

    # Precession of the RAAN
    # del_pRAAN will be in decimal degrees
    del_pRAAN = 360*(-3*J2*Re**2*n*t*cos(i)/(2*p0**2))
    pRAAN = RAAN + del_pRAAN

    # Precession of arg perigee
    # del_pAP will be in decimal degrees
    del_pAP = 360*(3*J2*Re**2*n*t*(5*cos(i)**2 - 1)/(4*p0**2))
    pAP = AP + del_pAP

    return [radians(pRAAN), radians(pAP)]


def calcArgLat(v, pAP):
    # Confirmed working
    # Inputs should be in radians, returns in radians
    Lat = pAP + v - 2*pi*(int((pAP + v)/(2*pi)))

    return Lat


def calcRADiff(sat, ArgLat):
    # Confirmed working
    # ArgLat should be in radians, Returns in radians
    i = radians(sat['INCLINATION'])

    if(0 < i < pi/2 and 0 < ArgLat < pi):
        RADiff = acos(cos(ArgLat)/sqrt(1-sin(i)**2 * sin(ArgLat)**2))
    elif(pi/2 < i < pi and pi < ArgLat < 2*pi):
        RADiff = acos(cos(ArgLat)/sqrt(1-sin(i)**2 * sin(ArgLat)**2))
    else:
        RADiff = 2*pi - acos(cos(ArgLat)/sqrt(1-sin(i)**2 * sin(ArgLat)**2))

    return RADiff


def calcGeocentricRADec(RADiff, pRAAN, ArgLat, manual_t=0):
    # Partially working, now accounts for rotation of earth
    # Inputs should be in radians, returns in radians
    GMST = utcDatetime2gmst(datetime.utcnow())*86400.0/86164.0
    earthrot = radians(GMST + manual_t*360.0)
    gcRA = -earthrot + RADiff + pRAAN - 2*pi*(int((RADiff + pRAAN)/(2*pi)))

    gcDec = (copysign(1, sin(ArgLat))) * acos(cos(ArgLat)/cos(RADiff))

    return [gcRA, gcDec]


def pol2cart(r, gcRA, gcDec):
    # Confirmed working
    # Takes inputs in radians, returns in km
    # Will give location of satellite

    Xg = r * cos(gcRA)*cos(gcDec)
    Yg = r * sin(gcRA)*cos(gcDec)
    Zg = r * sin(gcDec)

    return [Xg, Yg, Zg]


def LLA2cart(recvLat, recvLong, h):
    # Confirmed working
    # Lat long need to be in degrees, altitude in m
    # Returns in km

    a = 6378.137
    b = 6356.75231424518
    ecc = sqrt((a**2 - b**2)/a**2)
    N = a/sqrt(1 - ecc**2 * sin(recvLat)**2)
    X = (N + h) * cos(recvLat) * cos(recvLong)
    Y = (N + h) * cos(recvLat) * sin(recvLong)
    Z = (b**2/a**2 * N + h)*sin(recvLat)

    return [X, Y, Z]


def cart2RADec(satpos, obspos):
    # Returns in radians
    xg, yg, zg = satpos
    ag, bg, cg = obspos

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


def LLA2AzEl(satlat, satlong, obsvlat, obsvlong, satpos, obsvpos):

    satlat = radians(satlat)
    satlong = radians(satlong)

    dlong = satlong - obsvlong
    azimuth = atan2(sin(dlong)*cos(satlat), cos(obsvlat)*sin(satlat) - sin(obsvlat)*cos(satlat)*cos(dlong))
    azimuth = (azimuth + 2 * pi) % (2*pi)

    xd, yd, zd = satpos
    x, y, z = obsvpos

    dx = xd - x
    dy = yd - y
    dz = zd - z

    elevation = pi/2.0 - acos((x*dx + y*dy + z*dz) / sqrt((x**2+y**2+z**2)*(dx**2+dy**2+dz**2)))

    return [azimuth, elevation]


def ECEF2LLA(pos):
    # Confirmed working
    X = pos[0]
    Y = pos[1]
    Z = pos[2]

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
        if(i > 50):
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

    while(long > pi):
        long -= 2*pi
    while(long < -pi):
        long += 2*pi

    return [degrees(lati), degrees(long)]


def TLE2LLA(sat, obsvLLA, manual_t=0):
    t = epochDiff(sat) + manual_t/86400.0
    M = calcMA(sat, t)
    a = calcSMA(sat)
    v = calcTA(sat, M)
    P = calcPerigee(sat, a)
    pRAAN, pAP = precession(sat, a, t)
    ArgLat = calcArgLat(v, pAP)
    RADiff = calcRADiff(sat, ArgLat)
    r = geoDist(sat, P, v)
    gcRA, gcDec = calcGeocentricRADec(RADiff, pRAAN, ArgLat)
    cartcoords = pol2cart(r, gcRA, gcDec)
    obsvlat, obsvlong, height = obsvLLA

    obsvcoords = LLA2cart(obsvlat, obsvlong, height)

    alpha, delta, rg = cart2RADec(cartcoords, obsvcoords)
    LLAcoords = ECEF2LLA(cartcoords)
    satlat, satlong = LLAcoords
    azimuth, elevation = LLA2AzEl(satlat, satlong, obsvlat, obsvlong, cartcoords, obsvcoords)

    return [LLAcoords, azimuth, elevation, cartcoords]


def dataOutput(plot3d, plotRA, plotLLA, plotAz, duration, obsvLLA):

    plotting = max(plot3d, plotRA, plotLLA, plotAz)

    if(1 == plotting):
        path_to_chromedriver = 'C:\Users\Anthony\OneDrive\Programming\Sats\Sats\chromedriver.exe'
        browser = webdriver.Chrome(executable_path=path_to_chromedriver)
        url = 'http://www.n2yo.com/?s=25544'
        browser.get(url)
        sleep(10)
        tx = []
        ty = []
        tz = []
        RAfile = []
        latfile = []
        longfile = []
        azfile = []
        elfile = []
        i = 0

        t = epochDiff(sat)

        for time in xrange(1, duration, 5):
            output = TLE2LLA(sat, obsvLLA)
            satlat, satlong = output[0]
            az = degrees(output[1])
            el = degrees(output[2])
            cartcoords = output[3]

            # Lets webscrape
            n2lat = browser.find_element_by_id("satlat").text
            n2long = browser.find_element_by_id("satlng").text
            n2azimuth = browser.find_element_by_id("sataz").text
            n2el = browser.find_element_by_id("satel").text

            tx.append(cartcoords[0])
            ty.append(cartcoords[1])
            tz.append(cartcoords[2])
            latfile.append(satlat)
            longfile.append(satlong)
            azfile.append(abs(float(n2azimuth) - az))
            elfile.append(abs(float(n2el) - el))
            print i
            i += 5
            sleep(5)

    # Plots

    if(plot3d == 1):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        p = ax.plot(tx, ty, tz, color='black')

        # Draw sphere
        u = np.linspace(0, np.pi, 60)
        v = np.linspace(0, 2 * np.pi, 60)
        x = np.outer(np.sin(u), np.sin(v))
        y = np.outer(np.sin(u), np.cos(v))
        z = np.outer(np.cos(u), np.ones_like(v))
        ax.plot_wireframe(x*6371, y*6371, z*6371, color='blue')
        plt.show()

    if(plotRA == 1):
        plt.plot(RAfile)
        plt.ylabel('RA in rad')
        plt.xlabel('Time (s)')
        plt.show()

    if(plotLLA == 1):
        plt.plot(latfile, longfile)
        plt.ylabel('Lat (r)/ Long (b)')
        plt.xlabel('Time (s)')
        plt.show()

    if(plotAz == 1):
        t = []
        for i in range(0, len(azfile)):
            t.append(i)

        plt.plot(t, azfile, 'b', t, elfile, 'r')
        plt.ylabel('Azimuth error (deg)')
        plt.xlabel('Time (s)')
        plt.show()

#TLE = loadTLE()
#sat = TLE['ISS (ZARYA)']

#xc = []
#yc = []
#zc = []
"""for key in TLE:
    t1,t2,t3 = TLE2LLA(TLE[key])[1]
    xc.append(t1)
    yc.append(t2)
    zc.append(t3)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
p = ax.scatter(xc, yc, zc, color='black')

# Draw sphere
u = np.linspace(0, 2 * np.pi, 60)
v = np.linspace(0, 2 * np.pi, 60)
x = np.outer(np.sin(u), np.sin(v))
y = np.outer(np.sin(u), np.cos(v))
z = np.outer(np.cos(u), np.ones_like(v))
ax.plot_wireframe(x*6371, y*6371, z*6371,color='blue')
plt.show()
"""
#obsvLLA = [radians(-36.377518), radians(145.400044), 100]
#dataOutput(0, 0, 0, 0, obsvLLA, 6000)

#output = TLE2LLA(sat, obsvLLA)
#satlat, satlong = output[0]
#az = degrees(output[1])
#el = degrees(output[2])
