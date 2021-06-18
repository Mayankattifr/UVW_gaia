#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 31 23:52:13 2020

@author: mayank
"""

import numpy as np


def uvw(ra, dec, plx, pmra, pmde, rv):
    """
    :param ra: Right Ascension in degrees
    :param dec: Declination in degrees
    :param d: Distance in parsecs
    :param pmra: Proper motion in RA in milli-arcseconds/year
    :param pmde: Proper motion in Dec in milli-arcseconds/year
    :param rv: Radial velocity in km/s

    :return: U, V, W in km/s

    """

    k = 4.74047  # Equivalent of 1 A.U/yr in km/s
    A00 = -0.0548755604
    A01 = -0.8734370902
    A02 = -0.4838350155
    A10 = 0.4941094279
    A11 = -0.4448296300
    A12 = 0.7469822445
    A20 = -0.8676661490
    A21 = -0.1980763734
    A22 = 0.4559837762

    # Set as arrays in case ra, dec, etc were lists
    ra = np.array(ra)
    dec = np.array(dec)
    plx = np.array(plx)
    rv = np.array(rv)
    pmra = np.array(pmra)
    pmde = np.array(pmde)

    radcon = 3.1415926/180.0  # radian conversion factor
    cosd = np.array(np.cos(dec * radcon))
    sind = np.array(np.sin( dec * radcon))
    cosa = np.array(np.cos( ra * radcon))
    sina = np.array(np.sin( ra * radcon))

    vec1 = rv
    vec2 = k * pmra/plx
    vec3 = k * pmde/plx

    u = (A00*cosa*cosd + A01*sina*cosd + A02*sind) * vec1 + \
        (-A00*sina + A01*cosa) * vec2 + \
        (-A00*cosa*sind - A01*sina*sind + A02*cosd) * vec3
    v = (A10*cosa*cosd + A11*sina*cosd + A12*sind) * vec1 + \
        (-A10*sina + A11*cosa) * vec2 + \
        (-A10*cosa*sind - A11*sina*sind + A12*cosd) * vec3
    w = (A20*cosa*cosd + A21*sina*cosd + A22*sind) * vec1 + \
        (-A20*sina + A21*cosa) * vec2 + \
        (-A20*cosa*sind - A21*sina*sind + A22*cosd) * vec3
    u = u  # Flipping U to be positive towards Galactic center
    #return u+11.1, (v+12.24), (w+7.25)
    return u, (v), (w)	

def std2(u,m):
    u=np.array(u)
    n=len(u)
    ui=(u-m)**2
    ss=sum(ui/n)
    ss1=ss**0.5
    return ss1

	


