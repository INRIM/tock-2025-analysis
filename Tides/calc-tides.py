import numpy as np

import tintervals as ti

import os

import matplotlib.pyplot as plt
plt.ion()
plt.close("all")


def converter(s):
    return ti.iso2mjd(s + 'Z')


def calc_shift(dir_pot, dir_disp):
    # GEophysicist magic factor
    # see Kilian Stahl email on 25/3/2025
    magic_factor=1.3
    c = 299_792_458

    file_pot = os.path.join( dir_pot, "groundwater.app_pygtide.csv")
    file_disp = os.path.join( dir_disp, "groundwater.app_pygtide.csv")


    t_pot, pot = np.genfromtxt(file_pot, usecols=(0,1), converters={0:converter}, delimiter=',', skip_header=1, unpack=True)
    t_disp, disp = np.genfromtxt(file_disp, usecols=(0,1), converters={0:converter}, delimiter=',', skip_header=1, unpack=True)

    # check same timetags
    if not np.array_equal(t_pot, t_disp):
        raise ValueError("Please provide files with the same timetags")

    adjusted_shift = (magic_factor*pot - disp*9.81/1000)/c**2
    return t_pot, adjusted_shift



t_inrim, shift_inrim = calc_shift(dir_pot = "INRIM-potential", dir_disp = "INRIM-displacement")
t_ptb, shift_ptb = calc_shift(dir_pot = "PTB-potential", dir_disp = "PTB-displacement")


shift_diff = shift_inrim - shift_ptb
t = t_inrim

plt.figure()
plt.plot(t, shift_diff)

out = np.column_stack((t, shift_diff))
np.savetxt("INRIM-PTB.dat", out)