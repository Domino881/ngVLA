import sys
import math as m
import numpy as np
import scipy.constants
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

nfreqs = 2000 # number of points at which to sample frequency for output plots etc.

h_over_k = scipy.constants.Planck / scipy.constants.Boltzmann # in s*K
speedoflight = scipy.constants.speed_of_light

aeff_dict = {
#                     D  etaF_a etaF_b   epsp      epss     Ap    As     freqRange   
    "SKA":         [15.0,  0.92, 0.04, 280.0e-6, 154.0e-6, 0.89, 0.98,           []],
    "MeerKAT":     [13.5,  0.80, 0.04, 480.0e-6, 265.0e-6, 0.89, 0.98, [0.58, 3.05]],
    "Effelsberg":  [100.0,    0,    0,        0,        0,    0,    0, [1.00, 2.00]],
    "ngVLA":       [18.0,  0.80, 0.00, 280.0e-6, 154.0e-6, 0.89, 0.98, [1.00, 50.0]]
}

def heaviRange(x, range):
    if( range == [] ):
        return 1
    low, high = range
    # Inclusive of low, exclusive of high 
    return np.heaviside(x-low, 1) - np.heaviside(x-high, 1)

def get_aeff(telescope,plot):

    if( not aeff_dict.find(telescope) ):
        print("No idea what type of telescope I'm supposed to be calculating for. FAIL.")
        sys.exit(-1)

    D, etaF_a, etaF_b, epsp, epss, Ap, As, freqRange = aeff_dict[telescope]
    etaF = lambda freqGhz: etaF_a - etaF_b*freqGhz

    wavelength = scipy.constants.nu2lambda

    freq = np.logspace(np.log10(1.0), np.log10(50.0), 500)

    # Aperture efficiency
    delta   = 2*(Ap*epsp*epsp + As*epss*epss)**0.5
    DeltaPh = lambda freqGHz: 2*m.pi*delta/wavelength(freqGHz)
    etaPh = lambda freqGHz: np.exp(-(DeltaPh(freqGHz))**2.0) 
    etaD  = lambda freqGHz: 1.0 - 20.0*(wavelength(freqGHz)/D)**(1.5)

    # Overall aperture efficiency
    etaA  = lambda freqGHz: etaF(freqGHz)*etaPh(freqGHz)*etaD(freqGHz)

    if plot == True:
        plt.figure()
        plt.grid(True)
        plt.title("Aperture efficiency - %s dish"%(telescope))
        plt.ylabel("Aperture efficiency")
        plt.xlabel("Frequency (GHz)")
        plt.semilogx(freq, etaA(freq), 'o')
        plt.show()

    # Physical Collecting Area
    Aphys = m.pi*D*D/4.0

    # Effective collecting area
    Aeff = lambda freqGHz: Aphys * etaA(freqGHz) * heaviRange(freqGHz, freqRange)

    return Aeff

def get_tsys(telescope, gal, pwv, zenith, plot):

    # Receiver & Spillover Temperature
    if telescope == "SKA":
        Trcv = (
            lambda freqGHz: (15.0 + 30 * (freqGHz - 0.75) ** 2) * heaviRange(freqGHz, 0.35, 0.95)
            + (7.5) * heaviRange(freqGHz, 0.95, 4.60)
            + (4.4 + 0.69 * freqGHz) * heaviRange(freqGHz, 4.60, 50.0)
        )
        # Assumed to be this for all Bands but:
        # (a) is frequency dependent;
        # (b) is zenith angle dependent - 3 K is thought to be appropriate
        #     for zenith < 45 deg;
        # (c) the frequency dependence would actually be such that this should
        #     actually be a bit worse for Band 1 as it is not an octave feed.
        Tspill = lambda freqGHz: 3.0 + freqGHz*0.0 

    elif telescope == "MeerKAT":
        Trcv = (
            lambda freqGHz: (11.0 - 4.5 * (freqGHz - 0.58)) * heaviRange(freqGHz, 0.58, 1.02)
            + (7.5 + 6.8 * (np.abs(freqGHz - 1.65)) ** 1.5) * heaviRange(1.02, 1.65)
            + (7.5) * heaviRange(1.65, 3.05)
        )
        Tspill = lambda freqGHz: 4.0 + freqGHz*0.0

    elif telescope == "Effelsberg":
        Trcv = lambda freqGHz: 21.0 + freqGHz*0.0
        Tspill = lambda freqGHz: 0.0 + freqGHz*0.0 # don't know

    if telescope == "ngVLA":
        Trcv = lambda freqGHz: 20.0 + freqGHz*0.0  # placeholder value for now
        Tspill = lambda freqGHz: 1.0 + freqGHz*0.0 # placeholder value for now

    # Sky Temperature
    ## Tgal
    ## At the minute can only do 10th, 50th and 90the percentile values for Tgal
    ## Need to add any line of sight
    ## For now this is fine as it allows a direct comparison with Robert's calculations
    if (gal == "low"):
        tgal_pc = 10
        T408 = 17.1
    elif (gal == "medium"):
        tgal_pc = 50
        T408 = 25.2
    elif (gal == "high"):
        tgal_pc = 90
        T408 = 54.8

    # an off-plane approximation, need to do this more generally
    Tgal = lambda freqGHz: T408*(0.408/(freqGHz))**(2.75) 

    Tcmb = 2.73

    freq_array = np.genfromtxt("SKA_Tatm.txt", usecols=0)
    if (pwv == "low"):
        pwv_mm = 5.0
        Tatm_array = np.genfromtxt("SKA_Tatm.txt", usecols=1)
    elif (pwv == "medium"):
        pwv_mm = 10.0
        Tatm_array = np.genfromtxt("SKA_Tatm.txt", usecols=2)
    elif (pwv == "high"):
        pwv_mm = 20.0
        Tatm_array = np.genfromtxt("SKA_Tatm.txt", usecols=2)
    Tatm = interp1d(freq_array, Tatm_array, kind='cubic')
    Tsky = lambda freqGHz: Tgal(freqGHz) + Tcmb + Tatm(freqGHz)

    ### Opacity
    if (pwv == "low"):
        tau_array = np.genfromtxt("SKA_tau.txt", usecols=1)
    elif (pwv == "medium"):
        tau_array = np.genfromtxt("SKA_tau.txt", usecols=2)
    elif (pwv == "high"):
        tau_array = np.genfromtxt("SKA_tau.txt", usecols=2)
    tau = interp1d(freq_array, tau_array, kind='cubic')

    # tau = 0.01    # just eye-balled a reasonable value until I have the full function
    Tx = (
        lambda freqGHz, temp: (
            ((h_over_k * freqGHz * 1.0e9) / temp)
            / (np.exp((h_over_k * freqGHz * 1.0e9) / temp) - 1.0)
        )
        * temp
        * np.exp(tau(freqGHz) / np.cos(zenith * m.pi / 180.0))
    )

    Tsys = lambda f: Tx(f,(Trcv(f)+Tspill(f)+Tsky(f)))

    if (telescope == "ngVLA"):
        f = np.logspace(np.log10(1.0),np.log10(50.0),nfreqs)

    if (plot == True):
        plt.grid(True)
        plt.semilogx(f,Trcv(f),label='Receiver Temp.')
        plt.semilogx(f,Tspill(f),label='Spillover Temp.')
        plt.semilogx(f,Tsky(f),label='Sky Temp. (Gal+CMB+Atm)')
        plt.semilogx(f,Tsys(f),label='Tsys')
        plt.title("Temperature contributions, %2dth percentile $T_{\mathrm{Gal}}$, PWV %.1f mm"%(tgal_pc,pwv_mm))
        plt.ylabel("Temperature (K)")
        plt.xlabel("Frequency (GHz)")
        plt.legend()
        plt.show()

    return Tsys, f
