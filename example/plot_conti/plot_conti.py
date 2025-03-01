import numpy as np
from matplotlib import pyplot as plt
import matplotlib as mpl
from matplotlib.colors import ListedColormap
top = mpl.colormaps['Reds_r'].resampled(128)
bottom = mpl.colormaps['Blues'].resampled(128)
newcolors = np.vstack((top(np.linspace(0, 1, 128)),
                       bottom(np.linspace(0, 1, 128))))
residual_cmp = ListedColormap(newcolors, name='RedsBlue')
from matplotlib.patches import Ellipse
import os
from scipy import ndimage
from astropy.coordinates import SkyCoord

import sys
sys.path.insert(0,'../../')
from X22_model.disk_model import *
from radmc.setup import *
from radmc3dPy import *

sys.path.insert(0,'../')
from conti_fit_x22_builtin.find_center import find_center


filename = '~/project_data/cb68_edisk/CB68_SBLB_continuum_robust_0.0.image.tt0.fits'
ra_center, dec_center, disk_pa = find_center(filename)
# , x_lim=[2750,3250], y_lim=[2750,3250]
cb68_image = DiskImage(
    fname=filename,
    ra_deg = ra_center,
    dec_deg = dec_center,
    distance_pc = 140,
    rms_Jy = 21e-6, # convert to Jy/beam
    disk_pa = disk_pa,
    img_size_au = 60,
    remove_background=True
)

cb68_conti  = cb68_image.img
beam_pa     = cb68_image.beam_pa
beam_maj_au = cb68_image.beam_maj_au
beam_min_au = cb68_image.beam_min_au
# beam_area   = cb68_image.beam_area
au_per_pix  = cb68_image.au_per_pix
disk_pa     = cb68_image.disk_pa
size_au     = cb68_image.img_size_au
beam_area   = pi/(4*np.log(2)) * beam_maj_au*beam_min_au / (au_per_pix**2)


def radmc_conti(amax=0.1):
    model = radmc3d_setup(silent=False)
    model.get_mastercontrol(filename=None,
                            comment=None,
                            incl_dust=1,
                            incl_lines=1,
                            nphot=500000,
                            nphot_scat=100000,
                            scattering_mode_max=2,
                            istar_sphere=1,
                            num_cpu=None,
                            modified_random_walk=1
                            )
    model.get_continuumlambda(filename=None,
                            comment=None,
                            lambda_micron=None,
                            append=False)
    model.get_diskcontrol(  d_to_g_ratio=0.01,
                            a_max=amax, 
                            Mass_of_star=0.14, 
                            Accretion_rate=8e-7,
                            Radius_of_disk=30,
                            NR=200,
                            NTheta=200,
                            NPhi=10)
    model.get_heatcontrol(heat='accretion')
    os.system(f'radmc3d image npix {cb68_conti.shape[0]} sizeau {size_au} posang {disk_pa} incl 73 lambda 1300 noline')

    im = image.readImage()
    model_image = im.imageJyppix[:,:,0] /(140**2)
    I = ndimage.rotate(model_image, -2*disk_pa+beam_pa,reshape=False)
    sigmas = np.array([beam_maj_au, beam_min_au])/au_per_pix/(2*np.sqrt(2*np.log(2)))
    I = ndimage.gaussian_filter(I, sigma=sigmas)
    # rotate to align with image
    I = ndimage.rotate(I, -beam_pa,reshape=False)
    # convert to flux density in Jy/beam
    model_image = I*beam_area
    return model_image
# print(np.max(model_image))
# sigma_levels = [10, 20, 40, 80, 160, 320, 640]
# contour_levels = [level * 30e-6 for level in sigma_levels]

a_list = [1, 0.1, 0.01, 0.001]
fig, ax = plt.subplots(1, len(a_list), sharey=True, figsize=(4*len(a_list), 4))
fig.subplots_adjust(left=0.07, right=0.95, top=0.9, bottom=0.1, wspace=0.1, hspace=0.05)

for idx, a in enumerate(a_list):
    
    model_image = radmc_conti(amax=a)

    plot = ax[idx].imshow(model_image, cmap='Oranges', origin='lower', vmin=0.00, vmax=0.005)
    ax[idx].set_title(r'$a_{max} =$'+f'{a} mm')
    colorbar = fig.colorbar(plot, ax=ax[idx], pad=0.005, aspect=50, shrink=0.85)
    beam = Ellipse((12, 12), width=cb68_image.beam_min_au/cb68_image.au_per_pix, height=cb68_image.beam_maj_au/cb68_image.au_per_pix,
                angle=69, edgecolor='black', facecolor=None, lw=1.5, fill=False)
    ax[idx].add_patch(beam)


plt.savefig('compare_dif_amax.pdf', transparent=True)
os.system('make cleanall')
