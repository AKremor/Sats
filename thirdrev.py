import requests
from math import floor, pi, radians, degrees, sin, cos, acos, sqrt, copysign
from datetime import datetime,timedelta


def loadTLE():
    # Confirmed working
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

    nt = n*t
    Mt = M0 + 360*(nt - int(nt) - int((M0 + 360*(nt - int(nt)))/360))
    
    # Convert to from deg to rad
    Mt = radians(Mt)
    return Mt
    
def calcTA(sat,M):
    # Confirmed working
    # TODO: increase accuracy of this, maybe higher order or iterative
    # First we calculate eccentric anomaly, third order
    # This returns v in radians
    e = sat['ECCENTRICITY']
    E = M + e*sin(M) + e**2 *sin(M)*cos(M) + 1.0/2 * e**3 * sin(M) * (3*cos(M)**2 - 1)
    
    v = acos((cos(E)-e)/(1-e*cos(E)))
    
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
    # Non working
    n = sat['MNMOTION']
    e = sat['ECCENTRICITY']
    i = radians(sat['INCLINATION'])
    RAAN = sat['RAAN']
    AP = sat['ARGPERIGEE']
    # Radius of earth in km
    Re = 6378.135
    a1 = a/Re
    # Second gravity harmonic, J2 Propagation
    J2 =  0.0010826267
    d1 = (3*J2*(Re**2) * (3*cos(i)**2 - 1))/(4*a1**2 * ((1-e**2)**(3.0/2)))
    a0 = -1*a1*(134*(d1**3)/81 + d1**2 + d1/3 - 1)
    p0 = a0*(1-e**2) 
    
    print a1
    print J2
    print d1
    print a0
    print p0
    # Precession of the RAAN
    #THIS NEEDS TOBE VERIFIED AGAINST STEP 7, P SQARED MAY NOT BE IN DENOM
    del_pRAAN = 360*(-3*J2*Re**2*n*t*cos(i)/(2*p0**2))
    pRAAN = RAAN + del_pRAAN
    
    
    # Precession of arg perigee
    del_pAP = 360*(3*J2*Re**2*n*t*(5*cos(i)**2 - 1)/(4*p0**2))
    pAP = AP + del_pAP
    
    return [pRAAN, pAP]
    

def calcArgLat(v, pAP):
    # Confirmed working
    # Inputs should be in radians, returns in radians
    Lat = pAP + v - 360*(int( (pAP + v)/360))
    return Lat

    
def RADiff(sat, ArgLat):
    # Confirmed working
    # Returns in radians
    i = radians(sat['INCLINATION'])
    if((0<i<pi/2 and 0<ArgLat<pi) or (pi/2<i<pi and pi<ArgLat<2*pi)):
        del_RA = acos(cos(ArgLat)/sqrt(1-sin(i)**2 * sin(ArgLat)**2))
    else:
        del_RA = 360 - acos(cos(ArgLat)/sqrt(1-sin(i)**2 * sin(ArgLat)**2))
    
    return del_RA
    
def GeocentricRA(del_RA, pRAAN):
    
    # Confirmed working
    # Inputs should be in degrees, returns in degrees
    gcRA = del_RA + pRAAN - 360*(int((del_RA + pRAAN)/360))
    
    return gcRA
    
    
def GeocentricDec(ArgLat, del_RA):
    # Inputs in radians, output in radians
    gcDec = (copysign(1, sin(ArgLat))) * acos(cos(ArgLat)/cos(del_RA))
    return gcDec
    

def pol2cart(r, gcRA, gcDec):
    # Confirmed working
    # Takes inputs in deg, returns in km
    gcRA = radians(gcRA)
    gcDec = radians(gcDec)

    Xg = r * cos(gcRA)*cos(gcDec)
    Yg = r * sin(gcRA)*cos(gcDec)
    Zg = r * sin(gcDec)
    
    return [Xg, Yg, Zg]


TLE = loadTLE()
sat = TLE['SO-50']
t = epochDiff(sat)
M = calcMA(sat, t)
a = calcSMA(sat)
v = calcTA(sat,M)
P = calcPerigee(sat, a)
pRAAN = precession(sat,a,t)[0]
pAP = precession(sat,a,t)[1]
ArgLat = calcArgLat(v, pAP)
del_RA = RADiff(sat, ArgLat)
r = geoDist(sat, P, v)
#print('epochDiff',t)
#print('calcMA',degrees(calcMA(sat,t)))
#print('calcTA',degrees(calcTA(sat,M)))
#print('calcSMA',calcSMA(sat))
#print('calcPerigee',calcPerigee(sat,a))
#print('geoDist',geoDist(sat, P, v))
print('precession', precession(sat,a,t))
#print('ArgLat', ArgLat)
#print('RADiff', RADiff(sat, ArgLat))
print('GeocentricRA', GeocentricRA(degrees(del_RA), degrees(pRAAN)))
print('GeocentricDec', GeocentricDec(radians(ArgLat), radians(del_RA)))
print('pol2cart',pol2cart(7791.564499,14.9188694,44.73125163))