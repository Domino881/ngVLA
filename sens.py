#!/usr/local/bin/python3

#
# Author: Evan Keane
# Date: 27/09/2017
# Description: Determines SKA sensitivity curves etc. for various 
#              user-defined input configurations
#
# TODO
# 1. tidy up, add opacity, zenith angle and temp transforms etc.
# 2. read in array configuration CSV files
# 3. read in Haslam map for specific and average gl,gb values and ranges
# 4. add a check that the zenith angle is possible for the given sky coords!
#

# Load some useful packages
import argparse
import sys
import math as m
import numpy as np
from scipy import constants
import matplotlib.pyplot as plt
from Sensitivity.functions import *

h_over_k = 4.8e-11 # Planck's constant divided by Boltzmann's constant in s*K
kB = 1380.0        # Boltzmann's constant in Jy*m^2*K^-1
nfreqs = 2000      # number of points at which to sample frequency for output plots etc.

# Parse command line arguments
# parser = argparse.ArgumentParser()
# parser.add_argument('-radius', type=float, dest='radius', help='choose distance from the array centre, in km, for chosen sub-array (default: entire array)', default=150.0)
# parser.add_argument('-glgb', nargs=2, type=float, dest='coord', help='enter specific Galactic coordinates to use (default: gl=180.0, gb=-90.0) NOT DONE YET', default=[180.0,-90.0])
# parser.add_argument('-gallos', dest='gal', help='choose either 10th, 50th or 90th percentile value for the galaxy contribution to the sky temperature (low/medium/high, default: low)', default='low')
# parser.add_argument('-pwv', dest='pwv', help='choose either 5mm, 10mm or 20mm for the PWV value for choosing (a) the zenith opacity, and (b) the atmospheric temperature contribution to the sky temperature (low/medium/high, default: low)', default="low")
# parser.add_argument('-tel', dest='tel', help='Which telescopes to include - options are low, mk, ska (meaning SKA1-Mid dishes only), mid (meaning mk+ska), ngvla (default: ngvla)', default="ngvla")
# parser.add_argument('-nelements', type=int, dest='nelements', help='choose the inner nelements elements (default: entire array)', default=197)
# parser.add_argument('-o', dest='output', help='choose the type of output - plot, file or both (default: plot)', default="plot")
# parser.add_argument('-zenith', type=float, dest='zenith', help='choose a zenith angle in degrees (default: 0.0)', default=0.0)
# parser.add_argument('--version', action='version', version='%(prog)s 0.0.1')
# args = parser.parse_args()

# Set values from command line inputs
# tel = args.tel
# gl = args.coord[0]
# gb = args.coord[1]
# gal = args.gal
# output = args.output
# plot = True
# if output != "plot":
#     plot = False
# zenith = args.zenith
# pwv = args.pwv
# radius = args.radius
# nelements = args.nelements
radius = 150.0
gl, gb = [180.0,-90.0]
gal = 'low'
pwv = "low"
tel = "ngvla"
nelements = 197
output = "plot"
zenith = 0.0
plot=True

# Get the effective collecting area & system temperature
#if ( tel == "mid" or tel == "Mid" or tel == "MID" or ):
#Aeff_SKA = get_aeff("SKA",plot)
if ( tel == "mk" or tel == "MK" or tel == "meeerkat" or tel == "MeerKAT" ):
    Aeff  = get_aeff("MeerKAT",plot)
    Tsys, f  = get_tsys("MeerKAT",gal,pwv,zenith,plot)
elif ( tel == "ngvla" or tel == "ngVLA" ):
    Aeff = get_aeff("ngVLA",plot)
    Tsys, f = get_tsys("ngVLA",gal,pwv,zenith,plot)
    f = np.logspace(np.log10(1.0),np.log10(50.0),200)          # ngVLA freq range
else:
    print("No idea what telescope I'm supposed to be calculating for. FAIL.")
    sys.exit(-1)

fig, ax = plt.subplots()
plt.grid()
plt.semilogx(f,Aeff(f)/Tsys(f)) #(Trcv(f)+Tspill(f)+Tsky(f)))
fig.suptitle("Gain - single Dish - LoS dependent")
ax.set_ylabel("Aeff/Tsys (m^2/K)")
ax.set_xlabel("Frequency (GHz)")

# in K/Jy (LoS independent but no indication of Tsys impact on performance)
#f = np.logspace(np.log10(0.35),np.log10(50),200)
fig, ax = plt.subplots()
plt.grid()
plt.semilogx(f,Aeff(f)/(2*kB))
fig.suptitle("Gain - single Dish - LoS independent")
ax.set_ylabel("Aeff/(2*kB) (K/Jy)")
ax.set_xlabel("Frequency (GHz)")
plt.show()

SANITY_CHECK=False
if( SANITY_CHECK ):
    # Gain - user-requested sub-array
    # need to read in the configs, for now just do whole array
    # 
    # Also compare some actually relevant telescopes, maybe a different flag for imaging- and NIP-relevant ones to show
    f = np.logspace(np.log10(0.35),np.log10(50),200)
    fig, ax = plt.subplots()
    plt.grid()
    # plt.loglog(f,133.0*(Aeff_SKA(f)/Tsys_SKA(f))+64.0*(Aeff_MK(f)/Tsys_MK(f)), label='my numbers')
    fig.suptitle("Gain - entire SKA1-Mid array (133 SKA1 + 64 MeerKAT)")
    ax.set_ylabel("Aeff/Tsys (m^2/K)")
    ax.set_xlabel("Frequency (GHz)")

    # Sanity check
    roberts_freq = np.genfromtxt(sens_data_dir+"roberts_numbers_mid",usecols=0)
    roberts_gain = np.genfromtxt(sens_data_dir+"roberts_numbers_mid",usecols=2)
    plt.loglog(roberts_freq,roberts_gain,label='Robert numbers')
    plt.legend()

    Nska = 133
    Nmk  = 64
    if ((radius < 150.0) or (nelements < 197)):

        # Read in the array configuration
        array = np.genfromtxt("Configuration/MID_dist_metres.txt",dtype=[('name','S6'),('xm','f8'),('ym','f8')],skip_header=0,usecols=(1,2,3)) # NB this is an array of tuples, not a 2-D array, due to carrying the name label for each dish
        dist = np.zeros(np.size(array))    # work out the distance from centre in km
        for i in range(0,np.size(array)):
            dist[i] = 0.001*np.sqrt(array[i][1]*array[i][1]+array[i][2]*array[i][2])
        array = array[np.argsort(dist)]
        dist = dist[np.argsort(dist)]
        # Choose the sub-array of interest
        if (radius < 150.0): # if radius specified re-work out nelements, i.e. radius flag over-rules nelements flag if both set
            subarray = array[np.where( dist < radius)]
            subdist = dist[np.where( dist < radius)]
            nelements = subarray.shape[0]
        Nska = Nmk = 0
        for i in range(0,nelements):
            print(subarray[i][0][0], "bla")
            if subarray[i][0][0] == 77:   # ascii code for 'M' is 77
                Nmk +=1
            elif subarray[i][0][0] == 83: # ascii code for 'S' is 83
                Nska +=1
        if (tel == "mk"):
            Nska = 0
        if (tel == "ska"):
            Nmk = 0
        print("Considering a radius of %.1f"%(subdist[nelements-1]))
        print("Considering %d SKA and %d MeerKAT dishes"%(Nska,Nmk))
        fig, ax = plt.subplots()
        plt.grid()
        fig.suptitle("Gain - subarray radius %.1f km - %d SKA1 + %d MeerKAT"%(subdist[nelements-1],Nska,Nmk))
        ax.set_ylabel("Aeff/Tsys (m^2/K)")
        ax.set_xlabel("Frequency (GHz)")
        
        fig, ax = plt.subplots()
        plt.grid()
        fig.suptitle("Gain - subarray radius %.1f km - %d SKA1 + %d MeerKAT"%(subdist[nelements-1],Nska,Nmk))
        ax.set_ylabel("Aeff/(2*kB) (K/Jy)")
        ax.set_xlabel("Frequency (GHz)")

plt.show()