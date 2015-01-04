"""
This file provides the Airy functions

Airy wave is based on:
EM 1110-2-1100 WATER WAVE MECHANICS (Part II) 30 April 2002

Calculations shall be based on "Transitional Water"

"""

import numpy as np
from scipy.optimize import newton


def wl(t, d, g=9810.0):
    #
    # wave length based on "Transitional Water"
    #
    wave_l_func = lambda l: g * t ** 2 / 2.0 / np.pi * np.tanh(2.0 * np.pi * d / l) - l
    #
    wl0 = newton(wave_l_func, 10000.0)
    #
    assert isinstance(wl0, float)
    return wl0


# print wl(t, d)


def wk(h, t, d, z=0.0, l=0.0, phase=0.0, g=9810.0, trunc='false', output_params='full'):
    #
    # h      = wave height
    # t      = wave period
    # d      = water depth
    # z      = elevation of point of interest with respect to water surface
    # (-ve) for any point below water surface
    # l      = wave length
    # Phase  = Phase angle (radians)
    # g      = gravity accelerations in consistent units
    # trunc  = parameter to ignore wave inundation
    # true :
    #
    # if l is not provided then calculate l
    #
    if l == 0.0:
        l = wl(t, d, g)
    #
    # keep original value of z
    zc = z
    #
    # if z is above surface and less than wave amplitude, revert to z=0.0
    #
    if z > 0.0:
        if z <= h / 2.0:
            zc = 0.0
            #
    #
    # calculate the wave kinematics
    #
    el = h / 2.0 * np.cos(phase)
    vx = h / 2.0 * g * t / l * np.cosh(2.0 * np.pi * (zc + d) / l) / np.cosh(2.0 * np.pi * d / l) * np.cos(phase)
    vz = h / 2.0 * g * t / l * np.sinh(2.0 * np.pi * (zc + d) / l) / np.cosh(2.0 * np.pi * d / l) * np.sin(phase)
    ax = g * np.pi * h / l * np.cosh(2.0 * np.pi * (zc + d) / l) / np.cosh(2.0 * np.pi * d / l) * np.sin(phase)
    az = g * np.pi * h / l * np.sinh(2.0 * np.pi * (zc + d) / l) / np.cosh(2.0 * np.pi * d / l) * np.cos(phase)
    #
    # if z is below mudline or above wave amplitude, set values to 0.0
    #
    if d + z < 0.0:
        #
        el, vx, vz, ax, az = 0.0, 0.0, 0.0, 0.0, 0.0
    #
    #
    if trunc.lower() == "false":
        #
        limit = h / 2.0
    else:
        limit = el
    #
    #
    if z > limit:
        #
        el, vx, vz, ax, az = 0.0, 0.0, 0.0, 0.0, 0.0
    #
    #
    if output_params.lower() == 'el':
        #
        output = el
        #
    elif output_params.lower() == 'vx':
        #
        output = vx
        #
    elif output_params.lower() == 'vz':
        #
        output = vz
        #
    elif output_params.lower() == 'ax':
        #
        output = ax
        #
    elif output_params.lower() == 'az':
        #
        output = az
        #
    else:
        #
        output = el, vx, vz, ax, az
        #
    #
    return output



