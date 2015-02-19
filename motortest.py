# Stepper motor test
from motorControl import motor
from time import sleep

# Make a new motor
motor1 = motor( 7)

angle = 1
while(angle < 360):
    motor1.angle(180)
    angle += 1
    sleep(1)
