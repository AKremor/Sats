import stepperControl
import sensors


# Make the stepper
azimuth_stepper = stepperControl.Stepper(11, 12, 0.0001, 27 * 200/360.0)
elevation_stepper = stepperControl.Stepper(15, 16, 0.0001, 27 * 200/360.0)
i2c = sensors.init_compass(0, 0x1d)

print("Enter in the format 'direction steps', such as 'w 100' and press enter")
print("Press q-enter to quit")

while True:
    r_input = raw_input()
    vals = r_input.split()
    inkey = ord(vals[0])
    steps = int(vals[1])

    if inkey == 113:
        exit()
    elif inkey == 119:
        # up
        for i in range(steps):
            elevation_stepper.manual_step(True)
    elif inkey == 97:
        # left
        for i in range(steps):
            azimuth_stepper.manual_step(False)
    elif inkey == 100:
        # right
        for i in range(steps):
            azimuth_stepper.manual_step(True)
    elif inkey == 115:
        # down
        for i in range(steps):
            elevation_stepper.manual_step(False)

    real_azimuth = sensors.read_compass(i2c)
    print("Compass reads {}".format(real_azimuth))
