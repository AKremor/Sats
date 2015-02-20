# Stepper motor test
from math import radians
import satTLE
from motorControl import motor
from time import sleep

# Make a new motor
# Initial orientation should be 0 deg
azMotor = motor(7)
elMotor = motor(8) # Check pin
obsvLLA = [radians(-36.377518), radians(145.400044), 100]
satname = 'ISS (ZARYA)'
TLE = satTLE.loadTLE()

while True:
    az, el = satTLE.getAzEl(TLE, satname, obsvLLA)
    # .rotate takes the current angle
    azMotor.rotate(az)
    elMotor.rotate(el)
    sleep(1)
    