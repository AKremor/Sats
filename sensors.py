from math import sin, cos, atan2, radians, degrees, asin, sqrt
import time
import pynmea2
from smbus import SMBus
import serial


def init_GPS(port, baud):
    serial_port = serial.Serial(port, baud, timeout=0,
                                bytesize=serial.SEVENBITS,
                                parity=serial.PARITY_EVEN,
                                stopbits=serial.STOPBITS_ONE)

    return serial_port


def read_GPS(ser, out_queue):
    while True:
        time.sleep(1)
        line = ser.readline()

        try:
            msg = pynmea2.parse(line)

            # print('rkt lat {} lon {}'.format(msg.latitude, msg.longitude))
            out_queue.put(msg)
        except pynmea2.nmea.ParseError:
            pass


def rocket_solution(rocket, ground):

    azimuth = bearing(rocket.longitude, rocket.latitude,
                      ground.longitude, ground.latitude)
    distance = haversine(rocket.longitude, rocket.latitude,
                         ground.longitude, ground.latitude)
    elev = elevation(rocket.altitude, ground.altitude, distance)
    return azimuth, elev


def haversine(lon1, lat1, lon2, lat2):
    """
    Calculate the great circle distance between two points
    on the earth (specified in decimal degrees) in metres
    """
    # convert decimal degrees to radians
    lon1, lat1, lon2, lat2 = [radians(x) for x in [lon1, lat1, lon2, lat2]]

    # haversine formula
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
    c = 2 * asin(sqrt(a))
    r = 6371  # Radius of earth in kilometers. Use 3956 for miles
    return 1000 * c * r


def bearing(lon1, lat1, lon2, lat2):
    # convert decimal degrees to radians
    lon1, lat1, lon2, lat2 = [radians(x) for x in [lon1, lat1, lon2, lat2]]

    dlon = lon2 - lon1

    y = sin(sin(dlon) * cos(lat2))
    x = cos(lat1) * sin(lat2) - sin(lat1) * cos(lat2) * cos(dlon)

    return degrees(atan2(y, x))


def elevation(alt1, alt2, distance):
    # all inputs in m

    dalt = alt2 - alt1
    elev = atan2(dalt, distance)

    return degrees(elev)


def twos_comp_combine(msb, lsb):
    """Reproduced from post by davek,
    http://forum.pololu.com/viewtopic.php?f=32&t=9370"""

    twos_comp = 256*msb + lsb
    if twos_comp >= 32768:
        return twos_comp - 65536
    else:
        return twos_comp


def init_compass(bus_number, address):
    i2c = SMBus(bus_number)
    device_id = 0b1001001

    if i2c.read_byte_data(address, 0x0f) != device_id:
        raise IOError

    return i2c


def read_compass(i2c):
    LSM = 0x1d
    CTRL_6 = 0x25
    CTRL_7 = 0x26
    MAG_X_LSB = 0x08
    MAG_X_MSB = 0x09
    MAG_Y_LSB = 0x0A
    MAG_Y_MSB = 0x0B

    i2c.write_byte_data(LSM, CTRL_6, 0b00100000)  # set +/- 4 gauss full scale
    i2c.write_byte_data(LSM, CTRL_7, 0x00)  # turn on magnetometer

    magx = twos_comp_combine(i2c.read_byte_data(LSM, MAG_X_MSB),
                             i2c.read_byte_data(LSM, MAG_X_LSB))
    magy = twos_comp_combine(i2c.read_byte_data(LSM, MAG_Y_MSB),
                             i2c.read_byte_data(LSM, MAG_Y_LSB))

    angle = degrees(atan2(magy, magx))

    if angle < 0:
        angle += 360
    return angle
