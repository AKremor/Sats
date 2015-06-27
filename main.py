#!/Python27/python
print 'Content-type: text/html'
print

__author__ = 'Anthony'
import satTLE
import json


# Testing
# make a new satellite
TLE = satTLE.loadTLE()
ISS = satTLE.Satellite('ISS', TLE['ISS (ZARYA)'])

print json.dumps(ISS.LLAcoordinates(0)[0])
