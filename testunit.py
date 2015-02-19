from math import radians
import satTLE


obsvLLA = [radians(-36.377518), radians(145.400044), 100]
satname = 'ISS (ZARYA)'
TLE = satTLE.loadTLE()
for i in range(100):
    az, el = satTLE.returnAzEl(TLE,satname, obsvLLA)
    print az, el

    
    