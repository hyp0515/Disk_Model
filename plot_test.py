# import problem_setup as p
# import numpy as np
# from mpl_toolkits.mplot3d import axes3d
from matplotlib import pyplot as plt
# from matplotlib import cm
from radmc3dPy.image import *
from radmc3dPy.analyze import *


makeImage(npix=500,incl=30.,phi=0.,wav=100.,sizeau=120)   # This calls radmc3d 
fig2  = plt.figure()
a=readImage()
plotImage(a,log=True,au=True,maxlog=5,cmap='hot')


# makeImage(npix=500,incl=60.,phi=0.,wav=90.,sizeau=120)   # This calls radmc3d 
# fig2  = plt.figure()
# a=readImage()
# plotImage(a,log=True,au=True,maxlog=5,cmap='hot')

# os.system("radmc3d sed incl 90")

# # Plotting it "by hand", the SED as seen at 1 pc distance
# for incldeg in [0, 30, 45, 60, 75, 90]:
#     os.system(f"radmc3d sed incl {incldeg}")
#     fig3  = plt.figure()
#     s     = readSpectrum()
#     lam   = s[:,0]
#     nu    = 1e4*cc/lam
#     fnu   = s[:,1]
#     nufnu = nu*fnu
#     plt.plot(1e-9*nu,1e20*fnu,label=f'{incldeg}deg')
#     plt.xscale('log')
#     plt.yscale('log')
#     # plt.axis([1e-1, 5e+4, 1e-1, 1e4])
#     plt.xlim((1e+1, 5e+6))
#     plt.ylim(bottom=1e-1)
#     plt.xlabel('$\\nu [GHz]$')
#     plt.ylabel('$ Flux Density \; [mJy]$')
#     # plt.ylabel('$ Flux Density \; [\mathrm{erg}\,\mathrm{cm}^{-2}\,\mathrm{s}^{-1}\,\mathrm{Hz}^{-1}]$')
#     plt.legend()
#     plt.savefig(f'SEDs/mstar_50/{incldeg}deg.png')


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


# s = readSpectrum('spectrum.out')

# plt.plot(s[:,0], s[:,1])

# plt.show()