import numpy as np
from matplotlib import pyplot as plt
from radmc3dPy.image import *
from radmc3dPy.analyze import *
from disk_model import *
from problem_setup import problem_setup
from scipy.optimize import curve_fit

###############################################################################
# Plotter
def plot_sed(plot_nu=True, GHz=True, mjy=True):

    s = readSpectrum('spectrum.out')
    lam = s[:, 0]
    fnu = s[:, 1]/(140*140)
    if mjy is True:
        fnu = 1e26*fnu
        plt.ylabel('$ Flux Density \; [mJy]$')
        plt.ylim((1e-2, 1e4))
    else:
        fnu = 1e23*fnu
        plt.ylabel('$ Flux Density \; [Jy]$')
        plt.ylim((1e-5, 1e1))

    if plot_nu is True:
        nu = (1e-2*cc)/(1e-6*lam)
        if GHz is True:
            nu = 1e-9*nu
            plt.xlabel('$\\nu [GHz]$')
        else:
            plt.xlabel('$\\nu [Hz]$')
        fig = plt.plot(nu, fnu)
        plt.yscale('log')
        plt.xscale('log')
        plt.xlim((5e1, 1e3))
    else:
        fig = plt.plot(lam, fnu)
        plt.xlabel('$\\lambda [mm]$')
        plt.yscale('log')
        plt.xscale('log')
        plt.xlim((1e2, 1e6))
    return 
###############################################################################
"""
Different maximum grain sizes
"""
amax_list = [0.1, 0.01, 0.001]
for amax in amax_list:
    problem_setup(a_max=amax, Mass_of_star=0.14*Msun, Accretion_rate=0.14e-5*Msun/yr, Radius_of_disk=30*au, pancake=False)
    os.system("radmc3d sed incl 70")
    plot_sed()
label = [str(a*10)+' mm' for a in amax_list]
plt.legend(label)
observed_Freq = [233.8, 233.8, 233.8, 246.7, 246.7, 246.7, 246.7, 246.7]
observed_Flux = [   56,    55,    59,    62,    60,    60,    61,    66]
plt.scatter(observed_Freq, observed_Flux, color='black')
plt.savefig('different_amax')
plt.close()
###############################################################################
"""
Different mass of star
"""
mstar_list = [0.10, 0.15, 0.20, 0.25, 0.30]
for mstar in mstar_list:
    problem_setup(a_max=0.01, Mass_of_star=mstar*Msun, Accretion_rate=mstar*1e-5*Msun/yr, Radius_of_disk=30*au, pancake=False)
    os.system("radmc3d sed incl 70")
    plot_sed()
label = [str(mstar)+r' $M_{\odot}$' for mstar in mstar_list]
plt.legend(label)
observed_Freq = [233.8, 233.8, 233.8, 246.7, 246.7, 246.7, 246.7, 246.7]
observed_Flux = [   56,    55,    59,    62,    60,    60,    61,    66]
plt.scatter(observed_Freq, observed_Flux, color='black')
plt.savefig('different_mstar')
plt.close()

###############################################################################
"""
Different inclination angle
"""
incl_list = [0, 15, 45, 60, 75, 90]
problem_setup(a_max=0.01, Mass_of_star=0.14*Msun, Accretion_rate=0.14*1e-5*Msun/yr, Radius_of_disk=30*au, pancake=False)
for incl in incl_list:
    os.system(f"radmc3d sed incl {incl}")
    plot_sed()
label = [str(i)+r'$^{\circ}$' for i in incl_list]
plt.legend(label)
observed_Freq = [233.8, 233.8, 233.8, 246.7, 246.7, 246.7, 246.7, 246.7]
observed_Flux = [   56,    55,    59,    62,    60,    60,    61,    66]
plt.scatter(observed_Freq, observed_Flux, color='black')
plt.savefig('different_incl')
plt.close()

