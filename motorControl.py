#import RPi.GPIO as GPIO
from math import floor
from time import sleep
import time

class motor:
    def __init__(self, pin):
        #GPIO.setmode(GPIO.BOARD)
        self.pin = pin
        #GPIO.set(self.pin, OUT)
        self.currAngle = 0
        self.deltaAngle = 0


    def rotate(self, angle):
        #print str(angle)
        # The current azimuth of the sat
        self.deltaAngle = angle - self.currAngle
        
        STEPS_PER_ANGLE = 4*400/360.0
        steps = self.deltaAngle * STEPS_PER_ANGLE
        if steps < 1:
            return

        # This will slightly overstep
        while(steps > 0):
            print 'Moved ' + str(int(time.time())) + '  ' + str(steps) + '  ' + str(self.currAngle)
            #GPIO.output(self.pin, True)
            sleep(0.001)
            #GPIO.output(self.pin, False)
            sleep(0.001)
            steps -= 1
        
        self.currAngle = angle
        return 