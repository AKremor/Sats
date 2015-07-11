#!/Python27/python
print 'Content-type: text/html'
print

__author__ = 'Anthony'
import satTLE
import json
import time

TLE = satTLE.loadTLE()
ISS = satTLE.Satellite('ISS', TLE['ISS (ZARYA)'])

# For website output only
#print json.dumps(ISS.LLAcoordinates(0)[0])


def continuous_output():
    # Testing
    # make a new satellite
    TLE = satTLE.loadTLE()
    ISS = satTLE.Satellite('ISS', TLE['ISS (ZARYA)'])

    shepp = satTLE.Observer(-37.97357, 145.01636, 120)

    while True:
        satLLA, satCoords = ISS.LLAcoordinates(0)
        print shepp.getAzEl(satLLA, satCoords)
        time.sleep(5)


continuous_output()
