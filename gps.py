
import time
import serial
import pynmea2
from math import sin, cos, atan2, radians, degrees, asin, sqrt
from Queue import Queue

def init_GPS(port, baud):
	serial_port = serial.Serial(port, baud, timeout=0, bytesize=serial.SEVENBITS, parity=serial.PARITY_EVEN, stopbits=serial.STOPBITS_ONE)

	return serial_port

def read_GPS(ser, out_queue):
	while True:
		time.sleep(1)
		line = ser.readline()
		#print('Attempt {}'.format(line))
		try:
			msg = pynmea2.parse(line)
			#print(msg)
			print('rkt lat {} lon {}'.format(msg.latitude, msg.longitude))
			#gnd = pynmea2.parse("$GPGGA,184353.07,3625.1407,S,14540.7389,E,1,04,2.6,100.00,M,-33.9,M,,0000*6D")
			gnd = pynmea2.GGA('GP', 'GGA', ('184353.07', '3625.1407', 'S', '14540.7389', 'E', '1', '04' ,'2.6', '131.00', 'M', '-33.9', 'M','', '0000'))
			print('gnd lat {} lon {}'.format(gnd.latitude, gnd.longitude))
			print(rocket_solution(msg, gnd))
		except pynmea2.nmea.ParseError:
			pass


def rocket_solution(rocket, ground):

	# Here we find the altitude and heading
	r = 6000
	y = r * (rocket.latitude - ground.latitude)
	x = r * (rocket.longitude - ground.longitude)
	z = rocket.altitude - ground.altitude - (x**2 + y**2) /(2*r)
	#z = 100 - 50 - (x**2 + y**2)/(2*r)
	azimuth = atan2(x,y)
	elev = atan2(z, sqrt(x**2 + y**2))

	#print(azimuth)

	azimuth = bearing(rocket.longitude, rocket.latitude, ground.longitude, ground.latitude)
	distance = haversine(rocket.longitude, rocket.latitude, ground.longitude, ground.latitude)
	print(distance)
	print(rocket.altitude)
	print(ground.altitude)
	elev = elevation(rocket.altitude, ground.altitude, distance)
	return azimuth, elev

def haversine(lon1, lat1, lon2, lat2):
    """
    Calculate the great circle distance between two points
    on the earth (specified in decimal degrees) in metres
    """
    # convert decimal degrees to radians
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])

    # haversine formula
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
    c = 2 * asin(sqrt(a))
    r = 6371 # Radius of earth in kilometers. Use 3956 for miles
    return 1000 * c * r

def bearing(lon1, lat1, lon2, lat2):
	# convert decimal degrees to radians
	lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])

	dlon = lon2 - lon1
	dlat = lat2 - lat1
	y = sin(sin(dlon) * cos(lat2))
	x = cos(lat1) * sin(lat2) - sin(lat1) * cos(lat2) * cos(dlon)
	bearing = atan2(y, x)

	return degrees(bearing)

def elevation(alt1, alt2, distance):
	# all inputs in m

	dalt = alt2 - alt1
	elev = atan2(dalt, distance)

	return degrees(elev)
