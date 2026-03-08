#!/usr/bin/env python

""" Tables.py - Make tables  of atmospheric properties   (Python 3)

Adapted by
    Richard J. Kwan, Lightsaber Computing
from original programs by
    Ralph L. Carmichael, Public Domain Aeronautical Software

Revision History
Date         Vers Person Statement of Changes
2004 Oct 04  1.0  RJK    Initial program
2017 Jun 04  1.1  RLC    All indents are with spaces; all prints in ( )
                         New version for Python 3 that does integer div with //
2019 Sep 29  1.2  RLC    Changed polar radius of earth and renamed tables
2022 Mar 20  1.3  RLC    output files are *.output
"""

import sys, math

version = "1.3 (2022 March 20)"
greeting = "Tables - A Python program to compute atmosphere tables"
author = "Ralph L. Carmichael, Public Domain Aeronautical Software"
modifier = "Richard J. Kwan, Lightsaber Computing"
farewell = "Four files added to your directory."
finalmess = "Normal termination of tables."

#   P H Y S I C A L   C O N S T A N T S

FT2METERS = 0.3048      # mult. ft. to get meters (exact)
KELVIN2RANKINE = 1.8    # mult deg K to get deg R
PSF2NSM = 47.880258     # mult lb/sq.ft to get sq.m
SCF2KCM = 515.379       # mult slugs/cu.ft to get kg/cu.m
TZERO   = 288.15        # sea-level temperature, kelvins
PZERO   = 101325.0      # sea-level pressure, N/sq.m
RHOZERO = 1.225         # sea-level density, kg/cu.m
AZERO   = 340.294       # speed of sound at S.L.  m/sec
BETAVISC = 1.458E-6     # viscosity constant
SUTH    = 110.4         # Sutherland's constant, kelvins

def Atmosphere(alt):
    """ Compute temperature, density, and pressure in standard atmosphere.
    Correct to 86 km.  Only approximate thereafter.
    Input:
    alt geometric altitude, km.
    Return: (sigma, delta, theta)
    sigma   density/sea-level standard density
    delta   pressure/sea-level standard pressure
    theta   temperature/sea-level std. temperature
    """

    REARTH = 6356.766   # polar radius of the Earth (km)
    GMR = 34.163195     # hydrostatic constant (K/km)
    NTAB = 8            # length of tables

    htab = [ 0.0,  11.0, 20.0, 32.0, 47.0,
    51.0, 71.0, 84.852 ]
    ttab = [ 288.15, 216.65, 216.65, 228.65, 270.65,
    270.65, 214.65, 186.946 ]
    ptab = [ 1.0, 2.2336110E-1, 5.4032950E-2, 8.5666784E-3, 1.0945601E-3,
    6.6063531E-4, 3.9046834E-5, 3.68501E-6 ]
    gtab = [ -6.5, 0.0, 1.0, 2.8, 0, -2.8, -2.0, 0.0 ]

    h = alt*REARTH/(alt+REARTH) # geometric to geopotential altitude

    i=0; j=len(htab)-1
    while (j > i+1):
        k = (i+j)//2      # this is floor division in Python 3
        if h < htab[k]:
            j = k
        else:
            i = k
    tgrad = gtab[i]     # temp. gradient of local layer
    tbase = ttab[i]     # base  temp. of local layer
    deltah=h-htab[i]        # height above local base
    tlocal=tbase+tgrad*deltah   # local temperature
    theta = tlocal/ttab[0]  # temperature ratio

    if 0.0 == tgrad:
        delta=ptab[i]*math.exp(-GMR*deltah/tbase)
    else:
        delta=ptab[i]*math.pow(tbase/tlocal, GMR/tgrad)
    sigma = delta/theta
    return ( sigma, delta, theta )

def SimpleAtmosphere(alt):
    """ Compute temperature, density, and pressure in simplified
    standard atmosphere.

    Correct to 20 km.  Only approximate thereafter.

    Input:
        alt geometric altitude, km.
    Return: (sigma, delta, theta)
        sigma   density/sea-level standard density
        delta   pressure/sea-level standard pressure
        theta   temperature/sea-level std. temperature
    """
    REARTH = 6356.766   # polar radius of the Earth (km)
    GMR = 34.163195     # hydrostatic constant

    h = alt*REARTH/(alt+REARTH) # geometric to geopotential altitude

    if h<11.0:      # troposphere
        theta=(288.15-6.5*h)/288.15
        delta=math.pow(theta, GMR/6.5)
    else:       # stratosphere
        theta=216.65/288.15
        delta=0.2233611*math.exp(-GMR*(h-11.0)/216.65)
    sigma = delta/theta
    return ( sigma, delta, theta )

def MetricViscosity(theta):
    t=theta*TZERO
    return BETAVISC*math.sqrt(t*t*t)/(t+SUTH)
