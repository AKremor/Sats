from math import radians
from time import sleep
import satTLE


obsvLLA = [radians(-37.78661), radians(144.96849), 100]
satname = 'ISS (ZARYA)'
TLE = satTLE.loadTLE()
while True:
    az, el = satTLE.getAzEl(TLE,satname, obsvLLA)
    string = 'Azimuth: {:0.4f} Elevation: {:0.4f}'.format(az,el)
    print string
    sleep(1)

    
    