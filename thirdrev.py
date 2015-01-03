import requests
from math import floor, pi, radians, degrees, sin, cos, acos, sqrt, copysign,atan,asin,modf,atan2
from datetime import datetime,timedelta
from time import sleep
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

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

    UT = datetime2decHours(datetimeObj.time()) * 1.002737909
    T0 += UT

    GST = T0 % 24

    GST = GST * 15
    return GST


def datetime2decHours(time):
    # Confirmed working
	return time.hour + time.minute/60.0 + time.second/3600.0 + time.microsecond/3600000000.0


def loadTLE():
    # Confirmed working
    URL = 'http://www.amsat.org/amsat/ftp/keps/current/nasabare.txt'
    URL = 'http://www.celestrak.com/NORAD/elements/stations.txt'
    response = requests.get(URL)
    data = response.content
    TLE = {}

    # split the data into lines
    data = data.split('\n')
    
    for i in range(0, len(data)-2, 3):
        sat_name = data[i].rstrip()
        #sat_name = sat_name.replace(' ','_')
        #sat_name = sat_name.replace(')','')
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

    date_days = floor(float(epochDays))
    date_hours = floor((epochDays - date_days)*24)
    date_min = (epochDays - date_days - date_hours/24)*1440
    date_sec = (epochDays - date_days - date_hours/24 - date_min/1440)*86400
    date_epoch = str(epochYear) + ' ' + str(int(date_days)) + ' ' + str(int(date_hours)) + ' '+str(int(date_min))+ ' '+str(int(date_sec))
    epoch_time = datetime.strptime(date_epoch,"%Y %j %H %M %S")

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
    Mt = radians(Mt)
    return Mt


def calcTA(sat,M):
    # Confirmed working
    # TODO: increase accuracy of this, maybe higher order or iterative
    # First we calculate eccentric anomaly, third order
    # This returns v in radians
    e = sat['ECCENTRICITY']
    E = M + e*sin(M) + e**2 *sin(M)*cos(M) + 1.0/2 * e**3 * sin(M) * (3*cos(M)**2 - 1)

    x = sqrt(1-e)*cos(E/2)
    y = sqrt(1+e)*sin(E/2)
    v = 2*atan2(y,x)
    return v


def calcSMA(sat):
    # Confirmed working
    # Returns distance from earth center in km
    
    n = sat['MNMOTION']
    # mu is in km^3/day^2
    mu = 2.97554e15
    a = (mu/(2*pi*n)**2)**(1.0/3)
    return a


def calcPerigee(sat,a):
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


def precession(sat,a,t):
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
    J2 =  0.0010826267

    k2 = 1/2.0 * J2 * aE**2
    d1 = 3/2.0 * k2/a**2 *(3*cos(i)**2 - 1)/((1-e**2)**(1.0/3))

    a0 = a1*(1 - d1/3.0 - d1**2 - 134/81.0 * d1**3)
    p0 = a0 * (1-e**2)
    a0 = Re * a0
    p0 = Re * p0

    # Precession of the RAAN
    #del_pRAAN will be in decimal degrees
    del_pRAAN = 360*(-3*J2*Re**2*n*t*cos(i)/(2*p0**2))
    pRAAN = RAAN + del_pRAAN

    # Precession of arg perigee
    #del_pAP will be in decimal degrees
    del_pAP = 360*(3*J2*Re**2*n*t*(5*cos(i)**2 - 1)/(4*p0**2))
    pAP = AP + del_pAP

    return [radians(pRAAN), radians(pAP)]


def calcArgLat(v, pAP):
    # Confirmed working
    # Inputs should be in radians, returns in radians
    Lat = pAP + v - 2*pi*(int( (pAP + v)/(2*pi)))
    return Lat


def calcRADiff(sat, ArgLat):
    # Confirmed working
    # ArgLat should be in radians, Returns in radians
    i = radians(sat['INCLINATION'])

    if((0<i<pi/2 and 0<ArgLat<pi) or (pi/2<i<pi and pi<ArgLat<2*pi)):
        RADiff = acos(cos(ArgLat)/sqrt(1-sin(i)**2 * sin(ArgLat)**2))
    else:
        RADiff = 2*pi - acos(cos(ArgLat)/sqrt(1-sin(i)**2 * sin(ArgLat)**2))

    return RADiff


def calcGeocentricRA(RADiff, pRAAN):
    # Confirmed working
    # Inputs should be in radians, returns in radians
    gcRA = RADiff + pRAAN - 2*pi*(int((RADiff + pRAAN)/(2*pi)))

    return gcRA


def calcGeocentricDec(ArgLat, RADiff):
    # Inputs in radians, output in radians
    gcDec = (copysign(1, sin(ArgLat))) * acos(cos(ArgLat)/cos(RADiff))
    return gcDec


def pol2cart(r, gcRA, gcDec):
    # Confirmed working
    # Takes inputs in radians, returns in km
    # Will give location of satellite

    Xg = r * cos(gcRA)*cos(gcDec)
    Yg = r * sin(gcRA)*cos(gcDec)
    Zg = r * sin(gcDec)
    return [Xg, Yg, Zg]


def LLA2cart(recvLat=0, recvLong=0, recvAlt=100):
    # Confirmed working, note this doesn't agree with online converter because notECEF
    #lat long need tobe in degrees, altitude in m
    # Returns in km
    # Perform geodetic to geocentric conversions
    obsgcDec = radians(recvLat)
    
    # If not working, check whether the recvLong is needed 
    obsgcRA = utcDatetime2gmst(datetime.utcnow()) + recvLong
    
    while(obsgcRA > 360):
        obsgcRA -= 360
    obsgcRA = radians(obsgcRA)    

    # Calculate earth center to receiver location
    Re = 6378.135
    Rp = 6356.752

    rg = ((cos(obsgcDec)/Re)**2 + (sin(obsgcDec)/Rp)**2)**(-1.0/2) + recvAlt/1000
    #print('rg',rg)
    ag = rg * cos(obsgcRA)*cos(obsgcDec)
    bg = rg * sin(obsgcRA)*cos(obsgcDec)
    cg = rg * sin(obsgcDec)

    return [ag, bg, cg]


def cart2RADec(satpos, obspos):
    # Returns in radians
    xg = satpos[0]
    yg = satpos[1]
    zg = satpos[2]
    ag = obspos[0]
    bg = obspos[1]
    cg = obspos[2]

    xs = xg - ag
    ys = yg - bg
    zs = zg - cg

    if ys > 0 and xs > 0:
        alpha = atan(ys/xs)
    elif xs < 0:
        alpha = pi + atan(ys/xs)
    elif ys < 0 and xs > 0:
        alpha = 2*pi - atan(ys/xs)

    r = sqrt(xs**2 + ys**2 + zs**2)
    delta = asin(zs/r)
    #earthrot = radians(utcDatetime2gmst(datetime.utcnow()))
    #alpha = alpha + earthrot
    return [alpha,delta,r]

    
def RADec2AzAlt(alpha, delta, r, lat):   
    # alpha, delta in rad, lat in deg, r in km
    # Will return in radians
    HA = radians(utcDatetime2gmst(datetime.utcnow()))
    lat = radians(lat)
    
    Alt = asin(sin(delta)*sin(lat) + cos(delta)*cos(lat)*cos(HA))
    Az = acos((sin(delta) - sin(Alt)*sin(lat))/(cos(Alt)*cos(lat)))
    
    if sin(HA) < 0:
        Az = Az
    else:
        Az = 360 - Az
        
    return [Az, Alt]

def ECEF2LLA(pos,manual_t = 0):
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
    
    while(error > 0.0000001):
        if(i>50):
            print("enabled newton raphson modified guess")
            lati = atan(Z/sqrt(X**2 + Y**2))
            earthrot = radians(utcDatetime2gmst(datetime.utcnow()) + manual_t / 86400.0 * 360.0)
            #earthrot = 0
            long = long - earthrot
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

    earthrot = radians(utcDatetime2gmst(datetime.utcnow()) + manual_t / 86400.0 * 360.0)

    long = long - earthrot
    while(long > pi):
        long -= 2*pi
    while(long < -pi):
        long += 2*pi

    return [degrees(lati), degrees(long)]


# Plotting data

plot3d = 0
plotRA = 0
plotLLA = 0
plotting = max(plot3d, plotRA, plotLLA)

lat = -36.377518
TLE = loadTLE()
sat = TLE['ISS (ZARYA)']
t = epochDiff(sat)
#t = 1.7677141
MAfile = []
TAfile = []
tx = []
ty = []
tz = []
addflag = 0
RAfile = []
latfile = []
longfile= []
if(1==plotting):
    for time in range(1,12000):
        t += 1/86400.0
        M = calcMA(sat, t)
        a = calcSMA(sat)
        v = calcTA(sat,M)
        P = calcPerigee(sat, a)
        pRAAN = precession(sat,a,t)[0]
        pAP = precession(sat,a,t)[1]
        ArgLat = calcArgLat(v, pAP)
        RADiff = calcRADiff(sat, ArgLat)
        r = geoDist(sat, P, v)
        gcRA = calcGeocentricRA(RADiff, pRAAN)
        gcDec = calcGeocentricDec(ArgLat, RADiff)
        cartcoords = pol2cart(r, gcRA, gcDec)
        obsvcoords = LLA2cart(-36.377518, 145.400044, 100)
        RADec = cart2RADec(cartcoords, obsvcoords)
        alpha = RADec[0]
        delta = RADec[1]
        rg = RADec[2]
        AzAlt = RADec2AzAlt(alpha, delta, rg, lat)
        Az = AzAlt[0]
        Alt = AzAlt[1]
        LLAcoords = ECEF2LLA(cartcoords,t)
        print time
        
        
        
        if M<pi:
            addflag = 0

        MAfile.append(M + addflag*2*pi)
        TAfile.append(v)
        tx.append(cartcoords[0])
        ty.append(cartcoords[1])
        tz.append(cartcoords[2])
        RAfile.append(gcRA)
        latfile.append(LLAcoords[0])
        longfile.append(LLAcoords[1])


# Plots

timevalues = []
delta = []
for i in range(len(MAfile)):
    timevalues.append(i)
    delta.append(abs(MAfile[i]-TAfile[i]))

if(plot3d == 1):
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1,projection='3d')
    p = ax.plot(tx,ty,tz)
    plt.show()

if(plotRA == 1):
    plt.plot(RAfile)
    plt.ylabel('RA in rad')
    plt.xlabel('Time (s)')
    plt.show()

if(plotLLA == 1):
    plt.plot(timevalues, latfile, 'r', timevalues, longfile, 'b')
    #plt.plot(latfile,longfile)
    plt.ylabel('Lat (r)/ Long (b)')
    plt.xlabel('Time (s)')
    plt.show()


t = epochDiff(sat)
M = calcMA(sat, t)
a = calcSMA(sat)
v = calcTA(sat,M)
P = calcPerigee(sat, a)
pRAAN = precession(sat,a,t)[0]
pAP = precession(sat,a,t)[1]
ArgLat = calcArgLat(v, pAP)
RADiff = calcRADiff(sat, ArgLat)
r = geoDist(sat, P, v)
gcRA = calcGeocentricRA(RADiff, pRAAN)
gcDec = calcGeocentricDec(ArgLat, RADiff)
cartcoords = pol2cart(r, gcRA, gcDec)
obsvcoords = LLA2cart(-36.377518, 145.400044, 100)
RADec = cart2RADec(cartcoords, obsvcoords)
alpha = RADec[0]
delta = RADec[1]
rg = RADec[2]
AzAlt = RADec2AzAlt(alpha, delta, rg, lat)
Az = AzAlt[0]
Alt = AzAlt[1]
LLAcoords = ECEF2LLA(cartcoords)
#satpos = pol2cart(r, gcRA, gcDec)
#obspos = LLA2cart(-36.377518, 145.400044, 100)
print('epochDiff',t)
print('calcMA',degrees(calcMA(sat,t)))
print('calcTA',degrees(calcTA(sat,M)))
print('calcSMA',calcSMA(sat))
print('rg',r)
#print('calcPerigee',calcPerigee(sat,a))
print('pRAAN',pRAAN)
print('pAP',pAP)
print('geoDist',geoDist(sat, P, v))
print('precession', precession(sat,a,t))
print('ArgLat', ArgLat)
print('calcRADiff', calcRADiff(sat, ArgLat))
print('calcGeocentricRA',gcRA)
print('calcGeocentricDec',gcDec)
print('RA and Dec',alpha, delta)
print('pol2cart',cartcoords)
print('ECEF2LLA',LLAcoords)
#print('LLA2cart',obsvcoords)
#print('cart2RADec', degrees(alpha), degrees(delta),rg)
#print('RADec2AzAlt', degrees(Az),degrees(Alt))
#print('siderealtime',utcDatetime2gmst(datetime.utcnow()))
