import RPi.GPIO as GPIO
from math import floor, copysign
from time import sleep

class Stepper:
    def __init__(self, step_pin, dir_pin, sleep_pin = None):
        GPIO.setmode(GPIO.BOARD)
        self.step_pin = step_pin
        self.dir_pin = dir_pin
	self.sleep_pin = sleep_pin
        GPIO.setup(self.step_pin, GPIO.OUT)
        GPIO.setup(self.dir_pin, GPIO.OUT)
	if self.sleep_pin:
		GPIO.setup(self.sleep_pin, GPIO.OUT)
        self.current_angle = 0
	self.steps_to_make = 0

    def rotate_to_angle(self, move_to_angle, current_angle=None):

        steps_per_degree = 27 * 200/360.0

        if current_angle:
            # Angle from the compass to work with, zeroish error
            self.current_angle = current_angle

        # Calculate the angle we need to move through
        # Define CCW as +ve, CW -ve

        angle_delta = self.current_angle - move_to_angle
        steps_to_move = abs(angle_delta * steps_per_degree)

        if copysign(1, angle_delta) == 1:
            direction = 1 # FIX work out what direction this is
        else:
            direction = 0

        # Now perform the movement
	print 'Moving to angle {} from {} with {} steps+{} rem'.format(move_to_angle, self. current_angle, steps_to_move, self.steps_to_make)
        self.rotate_n_steps(steps_to_move, direction)

        # Assume rotation has worked, update current angle
        self.current_angle = move_to_angle

    def rotate_n_steps(self, num_steps, direction):
	
	self.steps_to_make += num_steps

        # Set the direction of rotation
        GPIO.output(self.dir_pin, direction)

        while 1 <= self.steps_to_make:
            GPIO.output(self.step_pin, True)
	    sleep(0.01)
            GPIO.output(self.step_pin, False)
            sleep(0.01)
	    self.steps_to_make -= 1

    def sleep(self, sleep_status):
	# Sleep is low, no sleep is high
	GPIO.output(self.sleep_pin, not sleep_status)
