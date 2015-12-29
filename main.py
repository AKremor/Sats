import threading
import gps
import satTLE
import time
import stepperControl
from Queue import Queue

TLE = satTLE.loadTLE()
ISS = satTLE.Satellite('ISS', TLE['EAGLE 2'])


def continuous_output():
    # Testing
    # make a new satellite
    TLE = satTLE.loadTLE()
    ISS = satTLE.Satellite('ISS', TLE['ISS (ZARYA)'])

    shepp = satTLE.Observer(-37.97357, 145.01636, 120)

    # Make the stepper
    azimuth_stepper = stepperControl.Stepper(11, 12)
    elevation_stepper = stepperControl.Stepper(7, 8)
    while True:
        satLLA, satCoords = ISS.LLAcoordinates(0)
        azimuth, elevation = shepp.getAzEl(satLLA, satCoords)

        # Now move the steppers
        azimuth_stepper.rotate_to_angle(azimuth)
	elevation_stepper.rotate_to_angle(elevation)
        time.sleep(5)


#continuous_output()

serial_port = gps.init_GPS('/dev/ttyAMA0', 38400)
q = Queue()

thread = threading.Thread(target=gps.read_GPS, args=(serial_port,q,))
thread.start()

