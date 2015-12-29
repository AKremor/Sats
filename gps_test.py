#!/usr/bin/python
import pynmea2
from serial import Serial
import time

serialPort = Serial("/dev/ttyAMA0", 38400, timeout=0)
if (serialPort.isOpen() == False):
	serialPort.open()
outStr = ''
inStr = ''
print('here1')
##### Set baud rate to 115200 #####
serialPort.flushInput()
serialPort.flushOutput()
print('here2')
outStr = '$PMTK314,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0*29\r\n'
 
while True:
	time.sleep(1)
	line = serialPort.readline()
	#print('Attempt {}'.format(line))
	try:
		msg = pynmea2.parse(line)
		print(msg)
	except pynmea2.nmea.ParseError:
		pass

