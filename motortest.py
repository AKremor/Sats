import RPi.GPIO as GPIO
GPIO.setmode(GPIO.BOARD)
GPIO.setup(11, GPIO.OUT)
GPIO.setup(12, GPIO.OUT)
GPIO.output(12, True)
from time import sleep


while True:
	GPIO.output(11, True)
	sleep(0.00010)
	print 'steppig'
	GPIO.output(11, False)
	sleep(00.00010)
