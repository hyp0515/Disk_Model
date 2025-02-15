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
from radmc3dPy import *


import sys
sys.path.insert(0,'../../')
from X22_model.disk_model import *
from radmc.setup import *
from CB68.data_dict import data_dict

sys.path.insert(0,'../')
from fit_with_GIdisk.find_center import find_center


# filename = '~/project_data/cb68_edisk/CB68_SBLB_continuum_robust_0.0.image.tt0.fits'
filename = data_dict["1.3_edisk"]["fname"]
ra_center, dec_center, disk_pa = find_center(filename)
# , x_lim=[2750,3250], y_lim=[2750,3250]
# ra_center  = '16:57:19.6428' # from 2d gaussian fitting
# dec_center = '-16:09:24.016'
# # Convert WCS to pixel coordinates
# coord_center = SkyCoord(ra=ra_center, dec=dec_center, unit=('hourangle', 'deg'))
# ra_center, dec_center = coord_center.ra.degree, coord_center.dec.degree

cb68_image = DiskImage(
    fname=filename,
    ra_deg = ra_center,
    dec_deg = dec_center,
    distance_pc = 140,
    rms_Jy = data_dict["1.3_edisk"]["sigma"], # convert to Jy/beam
    disk_pa = disk_pa,
    img_size_au = 120,
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
                        a_max=.1, 
                        Mass_of_star=0.14, 
                        Accretion_rate=4e-7,
                        Radius_of_disk=30,
                        NR=200,
                        NTheta=200,
                        NPhi=10)
model.get_heatcontrol(heat='accretion')
os.system(f'radmc3d image npix {cb68_conti.shape[0]} sizeau {size_au} incl 73 lambda 1300 noline')

im = image.readImage()
model_image = im.imageJyppix[:,:,0] /(140**2)
I = ndimage.rotate(model_image, -disk_pa+beam_pa,reshape=False)
sigmas = np.array([beam_maj_au, beam_min_au])/au_per_pix/(2*np.sqrt(2*np.log(2)))
I = ndimage.gaussian_filter(I, sigma=sigmas)
# rotate to align with image
I = ndimage.rotate(I, -beam_pa,reshape=False)
# convert to flux density in Jy/beam
model_image = I*beam_area
# print(np.max(model_image))
# sigma_levels = [10, 20, 40, 80, 160, 320, 640]
# contour_levels = [level * 30e-6 for level in sigma_levels]


fig, ax = plt.subplots(1, 3, sharey=True, figsize=(12, 4))
fig.subplots_adjust(left=0.07, right=0.95, top=0.9, bottom=0.1, wspace=0.1, hspace=0.05)

cb68 = ax[0].imshow(cb68_conti, cmap='Oranges', origin='lower', vmin=0.00)
# ax[0].plot([])
ax[0].set_title('CB68')
# ax[0].set_yticks([0, 49//2, 48])
# ax[0].set_yticklabels([-50, 0, 50])
# ax[0].set_ylabel('AU')
# ax[0].set_xticks([0, 49//2, 48])
# ax[0].set_xticklabels([-50, 0, 50])
# ax[0].set_xlabel('AU')
colorbar = fig.colorbar(cb68, ax=ax[0], pad=0.005, aspect=50, shrink=0.85)
beam = Ellipse((12, 12), width=cb68_image.beam_min_au/cb68_image.au_per_pix, height=cb68_image.beam_maj_au/cb68_image.au_per_pix,
               angle=69, edgecolor='black', facecolor=None, lw=1.5, fill=False)
ax[0].add_patch(beam)
# contours = ax[0].contour(cb68_conti, levels=contour_levels, colors='black', linewidths=0.8)
# colorbar.set_label('Intensity (Jy/beam)')

model = ax[1].imshow(model_image, cmap='Oranges', origin='lower', vmin=0.00)
ax[1].set_title('Synthetic Observation')
# ax[1].set_xticks([0, 49//2, 48])
# ax[1].set_xticklabels([-50, 0, 50])
# ax[1].set_xlabel('AU')
colorbar = fig.colorbar(model, ax=ax[1], pad=0.005, aspect=50, shrink=0.85)
beam = Ellipse((12, 12), width=cb68_image.beam_min_au/cb68_image.au_per_pix, height=cb68_image.beam_maj_au/cb68_image.au_per_pix,
               angle=69, edgecolor='black', facecolor=None, lw=1.5, fill=False)
ax[1].add_patch(beam)
# colorbar.set_label('Intensity (Jy/beam)')

residual = ax[2].imshow(cb68_conti-model_image, cmap=residual_cmp, origin='lower', vmin=-0.02, vmax=0.02)
ax[2].set_title('Residual')
# ax[2].set_xticks([0, 49//2, 48])
# ax[2].set_xticklabels([-50, 0, 50])
# ax[2].set_xlabel('AU')
colorbar = fig.colorbar(residual, ax=ax[2], pad=0.005, aspect=50, shrink=0.85)
colorbar.set_label('Intensity (Jy/beam)')

plt.savefig('residual.pdf', transparent=True)
os.system('make cleanall')
# print(np.unravel_index(np.argmax(convolved_image.T), convolved_image.shape))
# print(np.unravel_index(np.argmax(cb68_conti), cb68_conti.shape))

# cb68_alma_list = [
#     '/run/media/hyp0515/fd14f880-ba6f-450f-b82d-98ba3710dc5f/backup/CB68/cb68_alma/cb68_setup1_all.rob2.I.image.tt0.pbcor.smooth.dropdeg.fits',
#     '/run/media/hyp0515/fd14f880-ba6f-450f-b82d-98ba3710dc5f/backup/CB68/cb68_alma/cb68_setup2_all.rob2.I.image.tt0.pbcor.smooth.dropdeg.fits',
#     '/run/media/hyp0515/fd14f880-ba6f-450f-b82d-98ba3710dc5f/backup/CB68/cb68_alma/cb68_setup3_all.rob2.I.image.tt0.pbcor.smooth.dropdeg.fits'
# ]

# nu_list = [
#     233.8,
#     246.7,
#     95.0
# ]

# lam_list = cc * 1e4 / (np.array(nu_list) * 1e9)
# print(lam_list)

# sigma_list = [
#     30e-6,
#     40e-6,
#     26e-6
# ]

# observation_data = []
# beam_pa = []
# disk_posang = []
# beam_axis = []
# beam_area = []


# for i, fname in enumerate(cb68_alma_list):
#     ra_deg, dec_deg, disk_pa = find_center(fname)
#     image_class = DiskImage(
#         fname = fname,
#         ra_deg = ra_deg,
#         dec_deg = dec_deg,
#         distance_pc = 140,
#         rms_Jy = sigma_list[i], # convert to Jy/beam
#         disk_pa = disk_pa,
#         img_size_au = 100,
#         remove_background=True
#     )
#     observation_data.append(image_class.img)
#     beam_pa.append(image_class.beam_pa)
#     beam_axis.append([image_class.beam_maj_au, image_class.beam_min_au])
#     beam_area.append(image_class.beam_area)
#     disk_posang.append(disk_pa)

# size_au  = image_class.img_size_au
# au_per_pix  = image_class.au_per_pix
# npix = image_class.img.shape[0]
# pa = np.mean(disk_pa)


# with open('camera_wavelength_micron.inp', 'w+') as f:
#         f.write('%d\n'%(len(lam_list)))
#         for value in lam_list:
#             f.write('%13.6e\n'%(value))


# amax, Mstar, Mdot = -1, 0.14, np.log10(4e-7)



# model = radmc3d_setup(silent=False)
# model.get_mastercontrol(filename=None,
#                         comment=None,
#                         incl_dust=1,
#                         incl_lines=1,
#                         nphot=500000,
#                         nphot_scat=500000,
#                         scattering_mode_max=2,
#                         istar_sphere=1,
#                         num_cpu=1)
# model.get_linecontrol(filename=None,
#                     methanol='ch3oh leiden 0 0 0')
# model.get_continuumlambda(filename=None,
#                         comment=None,
#                         lambda_micron=None,
#                         append=False)
# model.get_diskcontrol(  d_to_g_ratio = 0.01,
#                         a_max=10**amax, 
#                         Mass_of_star=Mstar, 
#                         Accretion_rate=10**Mdot,
#                         Radius_of_disk=30,
#                         NR=200,
#                         NTheta=200,
#                         NPhi=10)
# model.get_heatcontrol(heat='accretion')


# os.system(f'radmc3d image npix {npix} sizeau {size_au} posang {pa} incl 73 loadlambda noline noscat')

# im = image.readImage()
# model_image_list = []
# for i in range(im.image.shape[2]):
#     model_image = im.image[:,:,i].T
#     I = ndimage.rotate(model_image, -disk_pa+beam_pa[i], reshape=False)
#     sigmas = np.array([beam_axis[i][0], beam_axis[i][1]])/au_per_pix/(2*np.sqrt(2*np.log(2)))
#     I = ndimage.gaussian_filter(I, sigma=sigmas)
#     # rotate to align with image
#     I = ndimage.rotate(I, -beam_pa[i], reshape=False)
#     # convert to flux density in Jy/beam
#     model_image = I*1e23*beam_area[i]
#     model_image_list.append(model_image)


    
# fig, ax = plt.subplots(1, 3)
# for i in range(len(model_image_list)):
#     print(np.max(model_image_list[i]))
#     ax[i].imshow(model_image_list[i])
# plt.show()
