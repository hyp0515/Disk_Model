import numpy as np
from matplotlib import pyplot as plt
import matplotlib as mpl
from matplotlib.colors import ListedColormap
top = mpl.colormaps['Reds_r'].resampled(128)
bottom = mpl.colormaps['Blues'].resampled(128)
newcolors = np.vstack((top(np.linspace(0, 1, 128)),
                       bottom(np.linspace(0, 1, 128))))
residual_cmp = ListedColormap(newcolors, name='RedsBlue')
import sys
sys.path.insert(0,'../../')
from X22_model.disk_model import *
import time
import astropy.table
import os

import astropy.constants as const
from astropy.io import fits
au = const.au.cgs.value
Msun = const.M_sun.cgs.value
yr = 365*24*3600
from astropy.coordinates import SkyCoord
from find_center import find_center
from matplotlib.patches import Ellipse

cb68_alma_list = [
    # '/run/media/hyp0515/fd14f880-ba6f-450f-b82d-98ba3710dc5f/backup/CB68/cb68_alma/cb68_setup1_all.rob2.I.image.tt0.pbcor.smooth.dropdeg.fits',
    # '/run/media/hyp0515/fd14f880-ba6f-450f-b82d-98ba3710dc5f/backup/CB68/cb68_alma/cb68_setup2_all.rob2.I.image.tt0.pbcor.smooth.dropdeg.fits',
    # '/run/media/hyp0515/fd14f880-ba6f-450f-b82d-98ba3710dc5f/backup/CB68/cb68_alma/cb68_setup3_all.rob2.I.image.tt0.pbcor.smooth.dropdeg.fits',
    '/run/media/hyp0515/fd14f880-ba6f-450f-b82d-98ba3710dc5f/backup/CB68_eDisk/CB68_SBLB_continuum_robust_0.0.image.tt0.fits'
]

lambda_list = [
    # 0.13,
    # 0.12,
    # 0.32,
    0.13,
]

vmax_list = [
    # 0.065,
    # 0.075,
    # 0.01,
    0.004
]

residual_list = [
    # 0.025,
    # 0.035,
    # 0.005,
    0.002
]

sigma_list = [
    # 30e-6,
    # 40e-6,
    # 26e-6,
    21e-6
]

# weight = [
#     0.475,
#     0.475,
#     0.05
# ]


# fig, axes = plt.subplots(4, figsize=(10, 10), sharex=True)

# a_array = np.linspace(-2, 1, 50, endpoint=True)

# labels    = ['log likelihood', r'$M_{*}$', r'$R_{d}$', r'$\dot{M}$']
# legends   = ['Setup 1', 'Setup 1+2', 'Setup 1+2+3']
# c = ['cyan', 'gold', 'green']
# for idx in range(len(c)):
#     f_alma_list = cb68_alma_list[:(idx+1)]
#     log_like  = []
#     msun_list = []
#     rd_list   = []
#     mdot_list = []
#     for amax in a_array:
#         opacity_table = generate_opacity_table_opt(
#             a_min=1e-6, a_max=10**amax, # min/max grain size
#             q=-3.5, # slope for dust size distribution, dn/da ~ a^q
#             dust_to_gas=0.01 # dust-to-gas ratio before sublimation
#         )
#         disk_property_table = generate_disk_property_table(opacity_table = opacity_table)
#         D = DiskFitting('CB68', opacity_table, disk_property_table)
#         cosI = np.deg2rad(73)
#         D.set_cosI(cosI=cosI)    

#         images = []
#         for i, fname in enumerate(f_alma_list):
#             ra_deg, dec_deg, disk_pa = find_center(fname)
#             DI_alma = DiskImage(
#                 fname = fname,
#                 ra_deg = ra_deg,
#                 dec_deg = dec_deg,
#                 distance_pc = 140,
#                 rms_Jy = sigma_list[i], # convert to Jy/beam
#                 disk_pa = disk_pa,
#                 img_size_au = 120,
#                 remove_background=True
#             )
#             D.add_observation(DI_alma, lambda_list[i])
#             images.append(DI_alma)    
#         D.fit()
#         log_likelihood = []
#         for i in range(len(f_alma_list)):
#             log_likelihood.append(-0.5 * np.sum((images[i].img - images[i].img_model) ** 2 / sigma_list[i]**2))
        
#         log_like.append(np.mean(log_likelihood))
#         msun_list.append(D.disk_model.Mstar/Msun)
#         rd_list.append(D.disk_model.Rd/au)
#         mdot_list.append(D.disk_model.Mdot/Msun*yr)
        

#     parms_list = [log_like, msun_list, rd_list, mdot_list]
    
#     for i in range(len(labels)):
#         ax = axes[i]
#         ax.plot(a_array, np.array(parms_list[i])or=c[idx], label=legends[idx])
#         if i == 0:
#             ax.legend(loc='lower right')
#         if i == 1:
#             ax.plot(a_array, 0.14*np.ones(a_array.shape), 'r:')
#         if i == 2:
#             ax.plot(a_array, 30*np.ones(a_array.shape), 'r:')
#         if i == 3:
#             ax.plot(a_array, 5e-7*np.ones(a_array.shape), 'r:')
#         ax.set_ylabel(labels[i])
        
# # plt.legend(loc='upper left')
# plt.savefig('log_likelihood.pdf', transparent=True)
# plt.show()





opacity_table = generate_opacity_table_opt(
    a_min=1e-6, a_max=.015, # min/max grain size
    q=-3.5, # slope for dust size distribution, dn/da ~ a^q
    dust_to_gas=0.01 # dust-to-gas ratio before sublimation
)
disk_property_table = generate_disk_property_table(opacity_table=opacity_table)
D = DiskFitting('CB68', opacity_table, disk_property_table)
cosI = np.cos(np.deg2rad(73))
D.set_cosI(cosI=cosI)    

images = []
pa = []
for i, fname in enumerate(cb68_alma_list):
    ra_deg, dec_deg, disk_pa = find_center(fname, x_lim=[2750,3250], y_lim=[2750,3250])
    ra_center  = '16:57:19.6428' # from 2d gaussian fitting
    dec_center = '-16:09:24.013'
    # Convert WCS to pixel coordinates
    coord_center = SkyCoord(ra=ra_center, dec=dec_center, unit=('hourangle', 'deg'))
    ra_deg, dec_deg = coord_center.ra.degree, coord_center.dec.degree
    DI_alma = DiskImage(
        fname = fname,
        ra_deg = ra_deg,
        dec_deg = dec_deg,
        distance_pc = 140,
        rms_Jy = sigma_list[i],
        img_size_au = 60,
        disk_pa=disk_pa,
        remove_background=True
    )
    D.add_observation(DI_alma, lambda_list[i])
    images.append(DI_alma)    
    pa.append(disk_pa)
D.fit(weights=None)
print('M_star [M_sun] =',D.disk_model.Mstar/Msun)
print('R_disk [au] =',D.disk_model.Rd/au)
print('M_dot [M_sun/yr] =', D.disk_model.Mdot/Msun*yr)


# fig, ax = plt.subplots(3, len(images), sharex=True, sharey=True, figsize=(15,5*len(images)))
# fig.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1, wspace=0.02, hspace=-0.36)

# for col in range(len(images)):
#     cb68 = ax[0, col].imshow(images[col].img, cmap='magma', origin='lower', vmin=-0.0005, vmax=vmax_list[col])
#     colorbar = fig.colorbar(cb68, ax=ax[0, col], pad=0.00, aspect=30, shrink=0.64)
#     if col == len(images)-1:
#         colorbar.set_label('Intensity (Jy/beam)')
        
#     ax[2, col].set_xticks([0, images[col].img.shape[0]//2, images[col].img.shape[0]-1])
#     ax[2, col].set_xticklabels([-100, 0, 100])
#     ax[2, col].set_xlabel('AU', fontsize=14)
#     # ax[0].contour(DI_alma.img, levels=[50*40e-6]ors='black', linewidths=0.65)
    
#     model = ax[1, col].imshow(images[col].img_model, cmap='magma', origin='lower', vmin=-0.0005, vmax=vmax_list[col])
#     colorbar = fig.colorbar(model, ax=ax[1, col], pad=0.00, aspect=30, shrink=0.64)
#     if col == len(images)-1:
#         colorbar.set_label('Intensity (Jy/beam)')
    
#     residual = ax[2, col].imshow(images[col].img-images[col].img_model, cmap=residual_cmp, origin='lower', vmin=-residual_list[col], vmax=residual_list[col])
#     colorbar = fig.colorbar(residual, ax=ax[2, col], pad=0.00, aspect=30, shrink=0.64)
#     if col == len(images)-1:
#         colorbar.set_label('Intensity (Jy/beam)')
    
# ax[0, 0].set_title('Setup 1 (1.3 mm)')
# try:
#     ax[0, 1].set_title('Setup 2 (1.2 mm)')
# except:
#     pass
# try:
#     # ax[0, 2].set_title('Setup 3 (3.2 mm)')
#     ax[0, 2].set_title('eDisk (1.3 mm)')
# except:
#     pass
# for row in range(len(images)):
#     ax[row, 0].set_yticks([])
    
# ax[0, 0].set_ylabel('Observation', fontsize=14)
# ax[1, 0].set_ylabel('Model', fontsize=14)
# ax[2, 0].set_ylabel('Residual', fontsize=14)
# fig, ax = plt.subplots(1,3, sharex=True, sharey=True, figsize=(15,5))
# fig.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1, wspace=0.0, hspace=0.0)
# for col in range(len(images)):
#     cb68 = ax[0].imshow(images[col].img, cmap='jet', origin='lower', vmin=-0.0005, vmax=vmax_list[col])
#     colorbar = fig.colorbar(cb68, ax=ax[0], pad=0.00, aspect=30, shrink=0.64)
#     if col == len(images)-1:
#         colorbar.set_label('Intensity (Jy/beam)')
        
#     beam = Ellipse((115, 10), width=DI_alma.beam_min_au/DI_alma.au_per_pix, height=DI_alma.beam_maj_au/DI_alma.au_per_pix,
#                angle=DI_alma.beam_pa, edgecolor='w', facecolor='w', lw=1.5, fill=True)
#     ax[0].add_patch(beam)
    
#     ax[2].set_xticks([0, images[col].img.shape[0]//2, images[col].img.shape[0]-1])
#     ax[2].set_xticklabels([-100, 0, 100])
#     ax[2].set_xlabel('AU', fontsize=14)
#     # ax[0].contour(DI_alma.img, levels=[50*40e-6]ors='black', linewidths=0.65)
    
#     model = ax[1].imshow(images[col].img_model, cmap='jet', origin='lower', vmin=-0.0005, vmax=vmax_list[col])
#     colorbar = fig.colorbar(model, ax=ax[1], pad=0.00, aspect=30, shrink=0.64)
#     if col == len(images)-1:
#         colorbar.set_label('Intensity (Jy/beam)')
    
#     residual = ax[2].imshow(images[col].img-images[col].img_model, cmap=residual_cmp, origin='lower', vmin=-residual_list[col], vmax=residual_list[col])
#     colorbar = fig.colorbar(residual, ax=ax[2], pad=0.00, aspect=30, shrink=0.64)
#     if col == len(images)-1:
#         colorbar.set_label('Intensity (Jy/beam)')
    
# ax[0].set_title('eDisk (1.3 mm)')
# for row in range(len(images)):
#     ax[row,].set_yticks([])
