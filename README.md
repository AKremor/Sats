Program will generate the azimuth and elevation for a given satellite from its TLE. This data is then used with an onboard compass and known GPS location to actuate an azimuth and elevation stepper to the appropriate pointing angles. From here the program refreshes the antenna position at 1Hz.

Also includes manual_motor.py which is a quick program to move the two steppers quickly.

Only main.py should need to be edited, specifically the Observer parameters to reflect your local coordinates.
Pinouts for the various components can be found in stepperControl.py and main.py
i2c needs to be enabled for the compass to work, if this is unavailable then the program can run relative to its starting angle by removing the real_azimuth parameter from the azimuth_stepper.
