"""Creates a stepper object, and controls the movement of the stepper."""
from time import sleep
import RPi.GPIO as GPIO


class Stepper(object):

    def __init__(self, step_pin, dir_pin, speed, steps_per_degree):
        GPIO.setmode(GPIO.BOARD)
        self.step_pin = step_pin
        self.dir_pin = dir_pin

        GPIO.setup(self.step_pin, GPIO.OUT)
        GPIO.setup(self.dir_pin, GPIO.OUT)

        self.current_angle = 0
        self.steps_to_make = 0
        self.speed = speed
        self.steps_per_degree = steps_per_degree
        self.step_threshold = 10

    def rotate_to_angle(self, move_to_angle, current_angle=None):

        if current_angle is not None:
            # Angle from the compass to work with, zeroish error
            self.current_angle = current_angle

        # Calculate the angle we need to move through
        # Define CCW as +ve, CW -ve

        angle_delta = self.current_angle - move_to_angle
        self.steps_to_make += angle_delta * self.steps_per_degree

        # Register the movement
        print 'To {}, from {}, with {} steps'.format(move_to_angle,
                                                     self.current_angle,
                                                     self.steps_to_make)
        self.current_angle -= angle_delta

    def has_steps(self):
        if (self.steps_to_make < self.step_threshold and
                self.steps_to_make > -self.step_threshold):
            return False

        return True

    def manual_step(self, direction):

        GPIO.output(self.dir_pin, direction)

        # Now make a step
        GPIO.output(self.step_pin, True)
        sleep(self.speed)
        GPIO.output(self.step_pin, False)
        sleep(self.speed)

    def step(self):

        if self.steps_to_make < 2 and self.steps_to_make > -2:
            return False

        if self.steps_to_make > 0:
            direction = -1
            GPIO.output(self.dir_pin, True)
        else:
            direction = 1
            GPIO.output(self.dir_pin, False)

        # Now make a step
        GPIO.output(self.step_pin, True)
        sleep(self.speed)
        GPIO.output(self.step_pin, False)
        sleep(self.speed)
        # print('step remaining {}'.format(self.steps_to_make))
        self.steps_to_make += direction

        return True
