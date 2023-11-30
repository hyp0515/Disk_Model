import numpy as np
from matplotlib import pyplot as plt
from radmc3dPy.image import *
from radmc3dPy.analyze import *
from disk_model import *
from vertical_profile_class import DiskModel_vertical
from problem_setup import problem_setup

from scipy.optimize import curve_fit


# problem_setup(a_max=1, Mass_of_star=0.5*Msun, Accretion_rate=5e-6*Msun/yr, Radius_of_disk=50*au)
# makeImage(npix=500,incl=45.,phi=0.,wav=100.,sizeau=120)   # This calls radmc3d 
# fig2  = plt.figure()
# a=readImage()
# plotImage(a,log=True,au=True,maxlog=5,cmap='hot')
# problem_setup(a_max=1, Mass_of_star=0.5*Msun, Accretion_rate=5e-6*Msun/yr, Radius_of_disk=50*au)
# makeImage(npix=500,incl=60.,phi=0.,wav=100.,sizeau=120)   # This calls radmc3d 
# fig2  = plt.figure()
# a=readImage()
# plotImage(a,log=True,au=True,maxlog=5,cmap='hot')
# problem_setup(a_max=1, Mass_of_star=.5*Msun, Accretion_rate=5e-6*Msun/yr, Radius_of_disk=50*au)
# for i in [105, 120, 135, 150, 165, 180]:    
#     makeImage(npix=500,incl=i,phi=0.,wav=80.,sizeau=120)
#     a=readImage()
#     plotImage(a,log=True,au=True,maxlog=5,cmap='hot')  


paperFreq = [6  , 14 , 43, 43 , 43 , 86, 99, 99, 99, 224, 340]
paperFlux = [0.7, 2.2, 10, 7.6, 2.3, 48, 58, 55, 4 , 256, 630]

# def sed(nu, k, alpha):
#     return k*nu**alpha

# params, _ = curve_fit(
#     sed,
#     paperFreq,
#     paperFlux
# )
# print(params)
plt.scatter(paperFreq, paperFlux, color='black')
# fre = np.linspace(0, 350, 100)
# plt.plot(fre, sed(fre, *params), color='blue')
# plt.show()

# # Plotting it "by hand", the SED as seen at 1 pc distance
# for amax_idx in [1, 0.1, 0.01]:
#     for msun_idx in [5, 1, 0.5, 0.1]:
#         msun = msun_idx*Msun
#         mdot = msun*1e-5/yr
#         problem_setup(a_max=amax_idx, Mass_of_star=msun, Accretion_rate=mdot, Radius_of_disk=50*au)
#         for incldeg in [0, 30, 45, 60, 75, 90]:
#             os.system(f"radmc3d sed incl {incldeg}")
#             fig3  = plt.figure()
#             s     = readSpectrum()
#             lam   = s[:,0]
#             nu    = 1e4*cc/lam
#             fnu   = s[:,1]
#             nufnu = nu*fnu
#             plt.plot(1e-9*nu,1e20*fnu,label=f'{incldeg}deg')
#             # plt.plot(paperFreq,paperFlux,'o')
#             plt.xscale('log')
#             plt.yscale('log')
#             # plt.axis([1e-1, 5e+4, 1e-1, 1e4])
#             plt.xlim((1e+1, 5e+6))
#             plt.ylim(bottom=1e-1)
#             plt.xlabel('$\\nu [GHz]$')
#             plt.ylabel('$ Flux Density \; [mJy]$')
#             # plt.ylabel('$ Flux Density \; [\mathrm{erg}\,\mathrm{cm}^{-2}\,\mathrm{s}^{-1}\,\mathrm{Hz}^{-1}]$')
#             plt.legend()
#             plt.savefig(f'SEDs/amax_{amax_idx}cm/mstar_{msun_idx}/{incldeg}deg.png')





# for amax_idx in [1, 0.1, 0.01, 0.001]:
#     problem_setup(a_max=amax_idx, Mass_of_star=1*Msun, Accretion_rate=1e-5*Msun/yr, Radius_of_disk=50*au)
#     os.system("radmc3d sed incl 0")
#     s     = readSpectrum()
#     lam   = s[:,0]
#     nu    = 1e4*cc/lam
#     fnu   = s[:,1]
#     nufnu = nu*fnu
#     plt.plot(1e-9*nu,1e20*fnu,label=f"{amax_idx}cm")
# plt.xscale('log')
# plt.yscale('log')
# plt.xlim((1e+1, 1e+7))
# plt.ylim(bottom=1e-1)
# # plt.ylim((1e-1,1e2))
# plt.xlabel('$\\nu [GHz]$')
# plt.ylabel('$ Flux Density \; [mJy]$')
# plt.legend()
# plt.title('Different maximum grain size')
# plt.savefig('different_amax')
# # plt.show()
# plt.close()

# for msun_idx in [20, 15, 10, 5, 1, 0.5]:
#     msun = msun_idx*Msun
#     mdot = msun*1e-5/yr
#     problem_setup(a_max=0.1, Mass_of_star=msun, Accretion_rate=mdot, Radius_of_disk=50*au)
#     os.system("radmc3d sed incl 75")
#     s     = readSpectrum()
#     lam   = s[:,0]
#     nu    = 1e4*cc/lam
#     fnu   = s[:,1]
#     nufnu = nu*fnu
#     plt.plot(1e-9*nu,1e20*fnu,label=f"{msun_idx}Msun")
# plt.xscale('log')
# plt.yscale('log')
# plt.xlim((1e+1, 5e+6))
# plt.ylim(bottom=1e-1)
# # plt.ylim((1e-1,1e3))
# plt.xlabel('$\\nu [GHz]$')
# plt.ylabel('$ Flux Density \; [mJy]$')
# plt.legend()
# plt.title('Different mass of star')
# plt.savefig('different_mstar')
# plt.close()

problem_setup(a_max=0.1, Mass_of_star=15*Msun, Accretion_rate=15e-5*Msun/yr, Radius_of_disk=50*au)
for incldeg in [0, 30, 45, 60, 75, 90]:
    os.system(f"radmc3d sed incl {incldeg}")
    s     = readSpectrum()
    lam   = s[:,0]
    nu    = 1e4*cc/lam
    fnu   = s[:,1]
    nufnu = nu*fnu
    plt.plot(1e-9*nu,1e20*fnu,label=f"{incldeg}deg")
plt.xscale('log')
plt.yscale('log')
plt.xlim((1e+1, 5e+6))
plt.ylim(bottom=1e-1)
plt.xlabel('$\\nu [GHz]$')
plt.ylabel('$ Flux Density \; [mJy]$')
plt.legend()
plt.title('Different inclined angle')
# plt.savefig('different_angle')
plt.show()
plt.close()





# problem_setup(a_max=0.1, Mass_of_star=0.1*Msun, Accretion_rate=0.1*1e-5*Msun/yr, Radius_of_disk=50*au)
# os.system("radmc3d sed incl 72")


# paperFreq = [6  , 14 , 43, 43 , 43 , 86, 99, 99, 99, 224, 340]
# paperFlux = [0.7, 2.2, 10, 7.6, 2.3, 48, 58, 55, 4 , 256, 630]
# plt.scatter(paperFreq,paperFlux,color='red')
# s     = readSpectrum()
# lam   = s[:,0]
# nu    = 1e4*cc/lam
# fnu   = s[:,1]
# nufnu = nu*fnu
# plt.plot(1e-9*nu,1e21*fnu)
# # plt.plot(paperFreq,paperFlux,'o')
# plt.xscale('log')
# plt.yscale('log')
# # plt.axis([1e-1, 5e+4, 1e-1, 1e4])
# plt.xlim((1e+0, 1e+3))
# plt.ylim(bottom=1e-1)
# plt.xlabel('$\\nu [GHz]$')
# plt.ylabel('$ Flux Density \; [mJy]$')
# # plt.ylabel('$ Flux Density \; [\mathrm{erg}\,\mathrm{cm}^{-2}\,\mathrm{s}^{-1}\,\mathrm{Hz}^{-1}]$')
# # plt.legend()
# # plt.savefig(f'SEDs/amax_{amax_idx}cm/mstar_{msun_idx}/{incldeg}deg.png')
# # plt.savefig('SED.png')
# plt.show()



# plt.xlim((3e+1, 3e+6))
# plt.xlabel('$\lambda\; [\mu \mathrm{m}$]')
# plt.ylabel('$\\nu F_\\nu \; [\mathrm{erg}\,\mathrm{cm}^{-2}\,\mathrm{s}^{-1}]$')

#
# Use the radmc3dPy.analyze tool set for plotting the SED, 
# this time let's plot nuLnu in units of Lsun
#
# fig4  = plt.figure()
# plotSpectrum(s,nulnu=True,lsun=True,xlg=True,ylg=False,micron=True)
# plt.axis([1e-1,1e4,1e-8, 50])
# plt.show()