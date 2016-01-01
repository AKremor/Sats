import threading
import sensors
import satTLE
import time
import stepperControl
from Queue import Queue


def rotate_antenna():
    serial_port = sensors.init_GPS('/dev/ttyAMA0', 38400)
    i2c = sensors.init_compass(0, 0x1d)
    q = Queue()

    t_gps = threading.Thread(target=sensors.read_GPS, args=(serial_port, q,))
    t_gps.start()

    sat_name = raw_input("Please enter the exact satellite name: ")
    # Testing
    # make a new satellite
    TLE = satTLE.loadTLE()
    ISS = satTLE.Satellite('ISS', TLE[sat_name])
    location = satTLE.Observer(-37.97357, 145.01636, 120)

    # Make the stepper
    azimuth_stepper = stepperControl.Stepper(11, 12, 0.001, 27 * 200/360.0)
    elevation_stepper = stepperControl.Stepper(15, 16, 0.001, 27 * 200/360.0)

    while True:

        satLLA, satCoords = ISS.LLAcoordinates(0)
        azimuth, elevation = location.getAzEl(satLLA, satCoords)

        real_azimuth = sensors.read_compass(i2c)
        # Register that we want to move
        azimuth_stepper.rotate_to_angle(azimuth, real_azimuth)
        elevation_stepper.rotate_to_angle(elevation)

        # Now step both motors
        if azimuth_stepper.has_steps():
            while azimuth_stepper.step():
                pass

        if elevation_stepper.has_steps():
            while elevation_stepper.step():
                pass

        time.sleep(1)


if __name__ == '__main__':
    rotate_antenna()
