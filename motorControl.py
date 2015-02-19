import RPi.GPIO as GPIO
from time import sleep

class motor:
    def __init__(self, pin):
        GPIO.setmode(GPIO.BOARD)
        self.pin = pin
        GPIO.set(self.pin, OUT)
        self.currAngle = 0
        self.deltaAngle = 0


    def rotate(angle):
        # The current azimuth of the sat
        self.deltaAngle = angle - self.currAngle
        
        STEPS_PER_ANGLE = 400/360.0
        steps = self.deltaAngle * STEPS_PER_ANGLE

        # This will slightly overstep
        while(steps > 0):
            GPIO.output(self.pin, True)
            sleep(0.01)
            GPIO.output(self.pin, False)
            steps -= 1