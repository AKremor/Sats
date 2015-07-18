#!/Python27/python
print 'Content-type: text/html'
print

__author__ = 'Anthony'
import satTLE
import json
import time
import stepperControl

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

    # Make the stepper
    azimuth_stepper = stepperControl.Stepper(7, 8)  # FIX

    while True:
        satLLA, satCoords = ISS.LLAcoordinates(0)
        azimuth, elevation = shepp.getAzEl(satLLA, satCoords)

        # Now move the stepper
        azimuth_stepper.rotate_to_angle(azimuth)
        time.sleep(5)


continuous_output()
