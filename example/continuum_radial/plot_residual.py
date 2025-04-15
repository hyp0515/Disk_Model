import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
import matplotlib.font_manager as fm
import matplotlib as mpl
from matplotlib.colors import ListedColormap
top = mpl.colormaps['Reds_r'].resampled(128)
bottom = mpl.colormaps['Blues'].resampled(128)

newcolors = np.vstack((top(np.linspace(0, 1, 128)),
                       bottom(np.linspace(0, 1, 128))))
residual_cmp = ListedColormap(newcolors, name='RedsBlue')
from scipy import ndimage
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d, interp2d
import warnings
import sys
import os
sys.path.insert(0,'../../')
from CB68.data_dict import data_dict

from radmc3dPy import image
from radmc3dPy.analyze import *
from make_conti import generate_model
sys.path.append("..")
from profile_radial.make_disk import generate_disk


from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs.utils import skycoord_to_pixel
from astropy.coordinates import SkyCoord
import astropy.constants as const
au = const.au.cgs.value
pc = const.pc.cgs.value


def crop_image(image, center, width_au, au_per_pix):
    width = int(width_au/au_per_pix)
    return image[center[0]-width//2:center[0]+width//2, center[0]-width//2:center[0]+width//2]

def rotate_image(image, posang):
    image = ndimage.rotate(image, posang, reshape=False, axes=(1, 0))
    image = np.nan_to_num(image, nan=0)
    return image

def radial_intensity(image_array, center, width):
  if center is None:
    peak_idx_x, peak_idx_y = np.unravel_index(np.argmax(image_array, axis=None), image_array.shape)
    center = peak_idx_y
  if width != 1:
    radial_profile = np.mean(image_array[:, center-width//2:center+width//2], axis=1)
  else:
    radial_profile = image_array[:, center]
  return radial_profile


def plot_residual(model):
    
    


    residual = edisk_image - model
    model_mask = model < 5*sigma
    mask = mask_cb68 | model_mask
    

    

    edisk_image[mask_cb68] = np.nan
    model[model_mask] = np.nan
    residual[mask] = np.nan
    
    chi_sq = np.nansum((residual)**2)

    fig, ax = plt.subplots(1,3, sharex=False, sharey=True, figsize=(15,5))
    fig.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.05, wspace=0.0, hspace=0.0)


    edisk = ax[0].imshow(edisk_image*1e3, origin='lower', cmap='magma', vmin=0, vmax=4.5)
    colorbar = fig.colorbar(edisk, ax=ax[0], pad=0.00, aspect=30, shrink=.87)
    ax[0].set_xticks([0, edisk_image.shape[0]//2, edisk_image.shape[0]-1])
    ax[0].set_xticklabels([-35, 0, 35])
    ax[0].set_xlabel('Offset [AU]', fontsize=14)
    ax[0].set_yticks([0, edisk_image.shape[0]//2, edisk_image.shape[0]-1])
    ax[0].set_yticklabels([-35, 0, 35])
    ax[0].set_ylabel('Offset [AU]', fontsize=14)
    ax[0].set_title('CB68', fontsize=16)

    model_ax = ax[1].imshow(model*1e3, origin='lower', cmap='magma', vmin=0, vmax=4.5)
    colorbar = fig.colorbar(model_ax, ax=ax[1], pad=0.00, aspect=30, shrink=.87)
    ax[1].set_xlabel('Offset [AU]', fontsize=14)
    ax[1].set_xticks([0, model.shape[0]//2, model.shape[0]-1])
    ax[1].set_xticklabels([-35, 0, 35], fontsize=14)
    ax[1].set_title('Model (Radiation)', fontsize=16)


    residual_ax = ax[2].imshow(residual*1e3, origin='lower', cmap = residual_cmp, vmin=-2, vmax=2)
    colorbar = fig.colorbar(residual_ax, ax=ax[2], pad=0.00, aspect=30, shrink=.87)
    colorbar.set_label(r'$I_{\nu}$ [mJy beam$^{-1}$]')
    ax[2].set_xlabel('Offset [AU]', fontsize=14)
    ax[2].set_xticks([0, residual.shape[0]//2, residual.shape[0]-1])
    ax[2].set_xticklabels([-35, 0, 35], fontsize=14)
    ax[2].set_title('Residual', fontsize=16)
    
    return chi_sq


distance_pc = 140
crop_sizeau = 80
sizeau = crop_sizeau
posang = 45

fits_data = data_dict["1.3_edisk"]
beam_axis = [fits_data["bmaj"], fits_data["bmin"]]
beam_area = beam_axis[0]*beam_axis[1]*np.pi/(4*np.log(2))
beam_pa = fits_data["bpa"]
sigma = fits_data["sigma"]

with fits.open(fits_data["fname"]) as hdul:
    data = hdul[0].data
    header = hdul[0].header
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        wcs = WCS(header=header)

au_per_pix_edisk = abs(header['CDELT1'])/180*np.pi*distance_pc*pc/au


edisk_image = crop_image(data, [2996, 2996], crop_sizeau, au_per_pix_edisk)
f_interp_edisk = interp2d(np.linspace(-(crop_sizeau//2), crop_sizeau//2, edisk_image.shape[0]),
                          np.linspace(-(crop_sizeau//2), crop_sizeau//2, edisk_image.shape[1]), edisk_image, kind='linear')
edisk_image = f_interp_edisk(np.linspace(-(crop_sizeau//2), crop_sizeau//2, 1000),
                            np.linspace(-(crop_sizeau//2), crop_sizeau//2, 1000))
peak_idx_x_edisk, peak_idx_y_edisk = np.unravel_index(np.argmax(edisk_image, axis=None), edisk_image.shape)
edisk_image = edisk_image[peak_idx_x_edisk-400:peak_idx_x_edisk+400, peak_idx_y_edisk-400:peak_idx_y_edisk+400]
mask_cb68 = edisk_image < 5*sigma



# a_list = [1e-2, 5e-2, 1e-1, 5e-1, 1e0, 1e1]
# L_star_list = [1e-1, 5e-1, 1e0, 3e0, 5e0, 1e1]
# Q_list = [0.1, 0.2, 0.3, 0.5, 1, 1.5]
# mdot_list = [1e-8, 1e-7, 1e-6, 1e-5]
# heat_list = ["radiation", "accretion"]

# model_im = image.readImage(fname=f'./simulation/outfile/conti_a_10.0_Lstar_0.1_Q_0.5_mdot_1e-06_accretion_scat.out')
# npix = model_im.nx
# pixel_area = (sizeau/npix/140)**2
# conv_image = model_im.imConv(dpc=140, fwhm=beam_axis, pa=-beam_pa)
# conv_image = conv_image.imageJyppix * beam_area/pixel_area/(140**2)
# f_interp_model = interp2d(np.linspace(-(crop_sizeau//2), crop_sizeau//2, conv_image.shape[0]),
#         np.linspace(-(crop_sizeau//2), crop_sizeau//2, conv_image.shape[1]), conv_image, kind='linear')
# conv_image = f_interp_model(np.linspace(-(crop_sizeau//2), crop_sizeau//2, 1000),
#                             np.linspace(-(crop_sizeau//2), crop_sizeau//2, 1000))
# peak_idx_x_model, peak_idx_y_model = np.unravel_index(np.argmax(conv_image, axis=None), conv_image.shape)
# conv_image = conv_image[peak_idx_x_model-400:peak_idx_x_model+400, peak_idx_y_model-400:peak_idx_y_model+400].T
# chisq = plot_residual(conv_image)

# plt.tight_layout()
# plt.savefig(f'./figures/residuals/accretion/a_10.0_Lstar_0.1_Q_0.5_mdot_1e-06_accretion.pdf', transparent=True)
# plt.close("all")




# model_im = image.readImage(fname=f'./simulation/outfile/conti_a_0.05_Lstar_5.0_Q_0.5_mdot_1e-06_radiation_scat.out')
# npix = model_im.nx
# pixel_area = (sizeau/npix/140)**2
# conv_image = model_im.imConv(dpc=140, fwhm=beam_axis, pa=-beam_pa)
# conv_image = conv_image.imageJyppix * beam_area/pixel_area/(140**2)
# f_interp_model = interp2d(np.linspace(-(crop_sizeau//2), crop_sizeau//2, conv_image.shape[0]),
#         np.linspace(-(crop_sizeau//2), crop_sizeau//2, conv_image.shape[1]), conv_image, kind='linear')
# conv_image = f_interp_model(np.linspace(-(crop_sizeau//2), crop_sizeau//2, 1000),
#                             np.linspace(-(crop_sizeau//2), crop_sizeau//2, 1000))
# peak_idx_x_model, peak_idx_y_model = np.unravel_index(np.argmax(conv_image, axis=None), conv_image.shape)
# conv_image = conv_image[peak_idx_x_model-400:peak_idx_x_model+400, peak_idx_y_model-400:peak_idx_y_model+400].T
# chisq = plot_residual(conv_image)

# plt.tight_layout()
# plt.savefig(f'./figures/residuals/radiation/a_0.05_Lstar_5.0_Q_0.5_mdot_1e-06_radiation.pdf', transparent=True)
# plt.close("all")




# a_list = [1e-2, 5e-2, 1e-1, 5e-1, 1e0, 1e1]
# L_star_list = [1e-1, 5e-1, 1e0, 3e0, 5e0, 1e1]
# Q_list = [0.1, 0.2, 0.3, 0.5, 1, 1.5]
# mdot_list = [1e-8, 1e-7, 1e-6]
# heat_list = ["radiation", "accretion"]



# chi_list_irr = []
# idx_list_irr = []

# chi_list_acc = []
# idx_list_acc = []



# for idx_a, a in enumerate(a_list):
#     for idx_l, L_star in enumerate(L_star_list):
#         for idx_q, Q in enumerate(Q_list):
#             for idx_mdot, mdot in enumerate(mdot_list):
#                 for heat in heat_list:
#                     try:
#                       model_im = image.readImage(fname=f'./simulation/outfile/conti_a_{a}_Lstar_{L_star}_Q_{Q}_mdot_{mdot}_{heat}_scat.out')
#                       npix = model_im.nx
#                       pixel_area = (sizeau/npix/140)**2
#                       conv_image = model_im.imConv(dpc=140, fwhm=beam_axis, pa=-beam_pa)
#                       conv_image = conv_image.imageJyppix * beam_area/pixel_area/(140**2)
#                       f_interp_model = interp2d(np.linspace(-(crop_sizeau//2), crop_sizeau//2, conv_image.shape[0]),
#                               np.linspace(-(crop_sizeau//2), crop_sizeau//2, conv_image.shape[1]), conv_image, kind='linear')
#                       conv_image = f_interp_model(np.linspace(-(crop_sizeau//2), crop_sizeau//2, 1000),
#                                                   np.linspace(-(crop_sizeau//2), crop_sizeau//2, 1000))
#                       peak_idx_x_model, peak_idx_y_model = np.unravel_index(np.argmax(conv_image, axis=None), conv_image.shape)
#                       conv_image = conv_image[peak_idx_x_model-400:peak_idx_x_model+400, peak_idx_y_model-400:peak_idx_y_model+400].T
#                       chisq = plot_residual(conv_image)

#                       plt.savefig(f'./figures/residuals/{heat}/a_{a}_Lstar_{L_star}_Q_{Q}_mdot_{mdot}_{heat}.pdf', transparent=True)
#                       plt.close("all")
#                       if heat == "radiation":
#                           chi_list_irr.append(chisq)
#                           idx_list_irr.append((idx_a, idx_l, idx_q, idx_mdot))
#                       else:
#                           chi_list_acc.append(chisq)
#                           idx_list_acc.append((idx_a, idx_l, idx_q, idx_mdot))
#                       print(chisq)
#                     except:
#                        pass
                    
# a_best_idx, L_star_best_idx, Q_best_idx, mdot_best_idx = idx_list_irr[chi_list_irr.index(min(chi_list_irr))]
# a_best = a_list[a_best_idx]
# L_star_best = L_star_list[L_star_best_idx]
# Q_best = Q_list[Q_best_idx]
# mdot_best = mdot_list[mdot_best_idx]
# print(f"Maximum grain size: {a_best} mm")
# print(f"Stellar luminosity: {L_star_best} Lsun")
# print(f"Toomre Q: {Q_best}")
# print(f"Mass accretion rate: {mdot_best} Msun/yr")


# a_best_idx, L_star_best_idx, Q_best_idx, mdot_best_idx = idx_list_acc[chi_list_acc.index(min(chi_list_acc))]
# a_best = a_list[a_best_idx]
# L_star_best = L_star_list[L_star_best_idx]
# Q_best = Q_list[Q_best_idx]
# mdot_best = mdot_list[mdot_best_idx]
# print(f"Maximum grain size: {a_best} mm")
# print(f"Stellar luminosity: {L_star_best} Lsun")
# print(f"Toomre Q: {Q_best}")
# print(f"Mass accretion rate: {mdot_best} Msun/yr")

# edisk_image = rotate_image(edisk_image, 45)


# plt.imshow(edisk_image, origin='lower')
# plt.contour(edisk_image, levels=[0.0001, 0.002, 0.003, 0.004, 0.005], colors='white')



a = 0.1
L_star = 5.0
Q = 0.5
mdot = 1e-6
heat = 'radiation'


model_im = image.readImage(fname=f'./simulation/outfile/conti_a_{a}_Lstar_{L_star}_Q_{Q}_mdot_{mdot}_{heat}_scat.out')
npix = model_im.nx
pixel_area = (sizeau/npix/140)**2
conv_image = model_im.imConv(dpc=140, fwhm=beam_axis, pa=-beam_pa)
conv_image = conv_image.imageJyppix * beam_area/pixel_area/(140**2)
f_interp_model = interp2d(np.linspace(-(crop_sizeau//2), crop_sizeau//2, conv_image.shape[0]),
        np.linspace(-(crop_sizeau//2), crop_sizeau//2, conv_image.shape[1]), conv_image, kind='linear')
conv_image = f_interp_model(np.linspace(-(crop_sizeau//2), crop_sizeau//2, 1000),
                            np.linspace(-(crop_sizeau//2), crop_sizeau//2, 1000))
peak_idx_x_model, peak_idx_y_model = np.unravel_index(np.argmax(conv_image, axis=None), conv_image.shape)
conv_image = conv_image[peak_idx_x_model-400:peak_idx_x_model+400, peak_idx_y_model-400:peak_idx_y_model+400].T
plot_residual(conv_image)

plt.savefig(f'./figures/residuals/{heat}/a_{a}_Lstar_{L_star}_Q_{Q}_mdot_{mdot}_{heat}.pdf', transparent=True)
plt.close("all")

# # model = generate_model(a, L_star, Q, mdot, heat = heat)
# model_im = image.readImage(fname=f'./simulation/outfile/conti_a_{a}_Lstar_{L_star}_Q_{Q}_mdot_{mdot}_{heat}_scat.out')
# npix = model_im.nx
# pixel_area = (sizeau/npix/140)**2
# conv_image = model_im.imConv(dpc=140, fwhm=beam_axis, pa=-beam_pa)
# conv_image = conv_image.imageJyppix * beam_area/pixel_area/(140**2)

# f_interp_model = interp2d(np.linspace(-(crop_sizeau//2), crop_sizeau//2, conv_image.shape[0]),
#                             np.linspace(-(crop_sizeau//2), crop_sizeau//2, conv_image.shape[1]), conv_image, kind='linear')
# conv_image = f_interp_model(np.linspace(-(crop_sizeau//2), crop_sizeau//2, 1000),
#                             np.linspace(-(crop_sizeau//2), crop_sizeau//2, 1000))
# peak_idx_x_model, peak_idx_y_model = np.unravel_index(np.argmax(conv_image, axis=None), conv_image.shape)
# conv_image = conv_image[peak_idx_x_model-400:peak_idx_x_model+400, peak_idx_y_model-400:peak_idx_y_model+400].T

# conv_image[mask] = np.nan

# residual = edisk_image - conv_image





# fig, ax = plt.subplots(1,3, sharex=False, sharey=True, figsize=(15,5))
# fig.subplots_adjust(left=0.05, right=0.97, top=0.9, bottom=0.1, wspace=0.0, hspace=0.0)


# edisk = ax[0].imshow(edisk_image*1e3, origin='lower', cmap='plasma', vmin=0, vmax=5)
# colorbar = fig.colorbar(edisk, ax=ax[0], pad=0.00, aspect=30, shrink=.98)
# # colorbar.set_label('Intensity (mJy/beam)')
# ax[0].set_xticks([0, edisk_image.shape[0]//2, edisk_image.shape[0]-1])
# ax[0].set_xticklabels([-np.round(sizeau//2), 0, np.round(sizeau//2)])
# ax[0].set_xlabel('Offset [AU]', fontsize=14)
# ax[0].set_yticks([0, edisk_image.shape[0]//2, edisk_image.shape[0]-1])
# ax[0].set_yticklabels([-np.round(sizeau//2), 0, np.round(sizeau//2)])
# ax[0].set_ylabel('Offset [AU]', fontsize=14)
# # ax[0].text(0.9, 0.9, 'eDisk (1.3 mm)', transform=ax[0].transAxes, fontsize=14, color='black')
# ax[0].set_title('CB68', fontsize=14)

# model = ax[1].imshow(conv_image*1e3, origin='lower', cmap='plasma', vmin=0, vmax=5)
# colorbar = fig.colorbar(model, ax=ax[1], pad=0.00, aspect=30, shrink=.98)
# # colorbar.set_label('Intensity (mJy/beam)')
# ax[1].set_xlabel('Offset [AU]', fontsize=14)
# ax[1].set_xticks([0, conv_image.shape[0]//2, conv_image.shape[0]-1])
# ax[1].set_xticklabels([-np.round(sizeau//2), 0, np.round(sizeau//2)])
# # ax[1].set_ylabel('Offset [AU]', fontsize=14)
# # ax[1].text(0.9, 0.9, 'eDisk (1.3 mm)', transform=ax[1].transAxes, fontsize=14, color='black')
# ax[1].set_title('Model (Radiation)', fontsize=14)


# residual_ax = ax[2].imshow(residual*1e3, origin='lower', cmap = residual_cmp, vmin=-2, vmax=2)
# colorbar = fig.colorbar(residual_ax, ax=ax[2], pad=0.00, aspect=30, shrink=.98)
# colorbar.set_label('Intensity (mJy/beam)')
# ax[2].set_xlabel('Offset [AU]', fontsize=14)
# ax[2].set_xticks([0, residual.shape[0]//2, residual.shape[0]-1])
# ax[2].set_xticklabels([-np.round(sizeau//2), 0, np.round(sizeau//2)])
# # ax[2].set_ylabel('Offset [AU]', fontsize=14)
# # ax[2].text(0.7, 0.7, 'Residual', transform=ax[2].transAxes, fontsize=14, color='black')
# ax[2].set_title('Residual', fontsize=14)
# plt.savefig('residual_radiation.pdf',transparent=True)
# plt.close("all")




# a = 0.1
# L_star = 3.0
# Q = 0.5
# mdot = 1e-5
# heat = 'radiation'


# # model = generate_model(a, L_star, Q, mdot, heat = heat)
# model_im = image.readImage(fname=f'./simulation/outfile/conti_a_{a}_Lstar_{L_star}_Q_{Q}_mdot_{mdot}_{heat}_scat.out')
# npix = model_im.nx
# pixel_area = (sizeau/npix/140)**2
# conv_image = model_im.imConv(dpc=140, fwhm=beam_axis, pa=-beam_pa)
# conv_image = conv_image.imageJyppix * beam_area/pixel_area/(140**2)

# f_interp_model = interp2d(np.linspace(-(crop_sizeau//2), crop_sizeau//2, conv_image.shape[0]),
#                             np.linspace(-(crop_sizeau//2), crop_sizeau//2, conv_image.shape[1]), conv_image, kind='linear')
# conv_image = f_interp_model(np.linspace(-(crop_sizeau//2), crop_sizeau//2, 1000),
#                             np.linspace(-(crop_sizeau//2), crop_sizeau//2, 1000))
# peak_idx_x_model, peak_idx_y_model = np.unravel_index(np.argmax(conv_image, axis=None), conv_image.shape)
# conv_image = conv_image[peak_idx_x_model-400:peak_idx_x_model+400, peak_idx_y_model-400:peak_idx_y_model+400].T

# conv_image[mask] = np.nan

# residual = edisk_image - conv_image





# fig, ax = plt.subplots(1,3, sharex=False, sharey=True, figsize=(15,5))
# fig.subplots_adjust(left=0.05, right=0.97, top=0.9, bottom=0.1, wspace=0.0, hspace=0.0)


# edisk = ax[0].imshow(edisk_image*1e3, origin='lower', cmap='plasma', vmin=0, vmax=5)
# colorbar = fig.colorbar(edisk, ax=ax[0], pad=0.00, aspect=30, shrink=.98)
# # colorbar.set_label('Intensity (mJy/beam)')
# ax[0].set_xticks([0, edisk_image.shape[0]//2, edisk_image.shape[0]-1])
# ax[0].set_xticklabels([-np.round(sizeau//2), 0, np.round(sizeau//2)])
# ax[0].set_xlabel('Offset [AU]', fontsize=14)
# ax[0].set_yticks([0, edisk_image.shape[0]//2, edisk_image.shape[0]-1])
# ax[0].set_yticklabels([-np.round(sizeau//2), 0, np.round(sizeau//2)])
# ax[0].set_ylabel('Offset [AU]', fontsize=14)
# # ax[0].text(0.9, 0.9, 'eDisk (1.3 mm)', transform=ax[0].transAxes, fontsize=14, color='black')
# ax[0].set_title('CB68', fontsize=14)

# model = ax[1].imshow(conv_image*1e3, origin='lower', cmap='plasma', vmin=0, vmax=5)
# colorbar = fig.colorbar(model, ax=ax[1], pad=0.00, aspect=30, shrink=.98)
# # colorbar.set_label('Intensity (mJy/beam)')
# ax[1].set_xlabel('Offset [AU]', fontsize=14)
# ax[1].set_xticks([0, conv_image.shape[0]//2, conv_image.shape[0]-1])
# ax[1].set_xticklabels([-np.round(sizeau//2), 0, np.round(sizeau//2)])
# # ax[1].set_ylabel('Offset [AU]', fontsize=14)
# # ax[1].text(0.9, 0.9, 'eDisk (1.3 mm)', transform=ax[1].transAxes, fontsize=14, color='black')
# ax[1].set_title('Model (Radiation)', fontsize=14)


# residual_ax = ax[2].imshow(residual*1e3, origin='lower', cmap = residual_cmp, vmin=-2, vmax=2)
# colorbar = fig.colorbar(residual_ax, ax=ax[2], pad=0.00, aspect=30, shrink=.98)
# colorbar.set_label('Intensity (mJy/beam)')
# ax[2].set_xlabel('Offset [AU]', fontsize=14)
# ax[2].set_xticks([0, residual.shape[0]//2, residual.shape[0]-1])
# ax[2].set_xticklabels([-np.round(sizeau//2), 0, np.round(sizeau//2)])
# # ax[2].set_ylabel('Offset [AU]', fontsize=14)
# # ax[2].text(0.7, 0.7, 'Residual', transform=ax[2].transAxes, fontsize=14, color='black')
# ax[2].set_title('Residual', fontsize=14)
# plt.savefig('residual_radiation_1e-5.pdf',transparent=True)
# plt.close("all")






# a = 1.0
# L_star = 1.0
# Q = 0.5
# mdot = 1e-6
# heat = 'accretion'


# # model = generate_model(a, L_star, Q, mdot, heat = heat)
# model_im = image.readImage(fname=f'./simulation/outfile/conti_a_{a}_Lstar_{L_star}_Q_{Q}_mdot_{mdot}_{heat}_scat.out')
# npix = model_im.nx
# pixel_area = (sizeau/npix/140)**2
# conv_image = model_im.imConv(dpc=140, fwhm=beam_axis, pa=-beam_pa)
# conv_image = conv_image.imageJyppix * beam_area/pixel_area/(140**2)

# f_interp_model = interp2d(np.linspace(-(crop_sizeau//2), crop_sizeau//2, conv_image.shape[0]),
#                             np.linspace(-(crop_sizeau//2), crop_sizeau//2, conv_image.shape[1]), conv_image, kind='linear')
# conv_image = f_interp_model(np.linspace(-(crop_sizeau//2), crop_sizeau//2, 1000),
#                             np.linspace(-(crop_sizeau//2), crop_sizeau//2, 1000))
# peak_idx_x_model, peak_idx_y_model = np.unravel_index(np.argmax(conv_image, axis=None), conv_image.shape)
# conv_image = conv_image[peak_idx_x_model-400:peak_idx_x_model+400, peak_idx_y_model-400:peak_idx_y_model+400].T
# mask_model = conv_image < 5*sigma

# residual = edisk_image - conv_image

# mask = mask_cb68 | mask_model

# edisk_image[mask_cb68] = np.nan
# conv_image[mask_model] = np.nan
# residual[mask] = np.nan





# fig, ax = plt.subplots(1,3, sharex=False, sharey=True, figsize=(15,5))
# fig.subplots_adjust(left=0.05, right=0.97, top=0.9, bottom=0.1, wspace=0.0, hspace=0.0)


# edisk = ax[0].imshow(edisk_image*1e3, origin='lower', cmap='plasma', vmin=0, vmax=5)
# colorbar = fig.colorbar(edisk, ax=ax[0], pad=0.00, aspect=30, shrink=.98)
# # colorbar.set_label('Intensity (mJy/beam)')
# ax[0].set_xticks([0, edisk_image.shape[0]//2, edisk_image.shape[0]-1])
# ax[0].set_xticklabels([-np.round(sizeau//2), 0, np.round(sizeau//2)])
# ax[0].set_xlabel('Offset [AU]', fontsize=14)
# ax[0].set_yticks([0, edisk_image.shape[0]//2, edisk_image.shape[0]-1])
# ax[0].set_yticklabels([-np.round(sizeau//2), 0, np.round(sizeau//2)])
# ax[0].set_ylabel('Offset [AU]', fontsize=14)
# # ax[0].text(0.9, 0.9, 'eDisk (1.3 mm)', transform=ax[0].transAxes, fontsize=14, color='black')
# ax[0].set_title('CB68', fontsize=14)

# model = ax[1].imshow(conv_image*1e3, origin='lower', cmap='plasma', vmin=0, vmax=5)
# colorbar = fig.colorbar(model, ax=ax[1], pad=0.00, aspect=30, shrink=.98)
# # colorbar.set_label('Intensity (mJy/beam)')
# ax[1].set_xlabel('Offset [AU]', fontsize=14)
# ax[1].set_xticks([0, conv_image.shape[0]//2, conv_image.shape[0]-1])
# ax[1].set_xticklabels([-np.round(sizeau//2), 0, np.round(sizeau//2)])
# # ax[1].set_ylabel('Offset [AU]', fontsize=14)
# # ax[1].text(0.9, 0.9, 'eDisk (1.3 mm)', transform=ax[1].transAxes, fontsize=14, color='black')
# ax[1].set_title('Model (Accretion)', fontsize=14)


# residual_ax = ax[2].imshow(residual*1e3, origin='lower', cmap = residual_cmp, vmin=-2, vmax=2)
# colorbar = fig.colorbar(residual_ax, ax=ax[2], pad=0.00, aspect=30, shrink=.98)
# colorbar.set_label('Intensity (mJy/beam)')
# ax[2].set_xlabel('Offset [AU]', fontsize=14)
# ax[2].set_xticks([0, residual.shape[0]//2, residual.shape[0]-1])
# ax[2].set_xticklabels([-np.round(sizeau//2), 0, np.round(sizeau//2)])
# # ax[2].set_ylabel('Offset [AU]', fontsize=14)
# # ax[2].text(0.7, 0.7, 'Residual', transform=ax[2].transAxes, fontsize=14, color='black')
# ax[2].set_title('Residual', fontsize=14)
# plt.savefig('residual_accretion.pdf',transparent=True)
# plt.close("all")





# a = 10.0
# L_star = 0.1
# Q = 0.5
# mdot = 1e-5
# heat = 'accretion'


# # model = generate_model(a, L_star, Q, mdot, heat = heat)
# model_im = image.readImage(fname=f'./simulation/outfile/conti_a_{a}_Lstar_{L_star}_Q_{Q}_mdot_{mdot}_{heat}_scat.out')
# npix = model_im.nx
# pixel_area = (sizeau/npix/140)**2
# conv_image = model_im.imConv(dpc=140, fwhm=beam_axis, pa=-beam_pa)
# conv_image = conv_image.imageJyppix * beam_area/pixel_area/(140**2)

# f_interp_model = interp2d(np.linspace(-(crop_sizeau//2), crop_sizeau//2, conv_image.shape[0]),
#                             np.linspace(-(crop_sizeau//2), crop_sizeau//2, conv_image.shape[1]), conv_image, kind='linear')
# conv_image = f_interp_model(np.linspace(-(crop_sizeau//2), crop_sizeau//2, 1000),
#                             np.linspace(-(crop_sizeau//2), crop_sizeau//2, 1000))
# peak_idx_x_model, peak_idx_y_model = np.unravel_index(np.argmax(conv_image, axis=None), conv_image.shape)
# conv_image = conv_image[peak_idx_x_model-400:peak_idx_x_model+400, peak_idx_y_model-400:peak_idx_y_model+400].T

# conv_image[mask] = np.nan

# residual = edisk_image - conv_image





# fig, ax = plt.subplots(1,3, sharex=False, sharey=True, figsize=(15,5))
# fig.subplots_adjust(left=0.05, right=0.97, top=0.9, bottom=0.1, wspace=0.0, hspace=0.0)


# edisk = ax[0].imshow(edisk_image*1e3, origin='lower', cmap='plasma', vmin=0, vmax=5)
# colorbar = fig.colorbar(edisk, ax=ax[0], pad=0.00, aspect=30, shrink=.98)
# # colorbar.set_label('Intensity (mJy/beam)')
# ax[0].set_xticks([0, edisk_image.shape[0]//2, edisk_image.shape[0]-1])
# ax[0].set_xticklabels([-np.round(sizeau//2), 0, np.round(sizeau//2)])
# ax[0].set_xlabel('Offset [AU]', fontsize=14)
# ax[0].set_yticks([0, edisk_image.shape[0]//2, edisk_image.shape[0]-1])
# ax[0].set_yticklabels([-np.round(sizeau//2), 0, np.round(sizeau//2)])
# ax[0].set_ylabel('Offset [AU]', fontsize=14)
# # ax[0].text(0.9, 0.9, 'eDisk (1.3 mm)', transform=ax[0].transAxes, fontsize=14, color='black')
# ax[0].set_title('CB68', fontsize=14)

# model = ax[1].imshow(conv_image*1e3, origin='lower', cmap='plasma', vmin=0, vmax=5)
# colorbar = fig.colorbar(model, ax=ax[1], pad=0.00, aspect=30, shrink=.98)
# # colorbar.set_label('Intensity (mJy/beam)')
# ax[1].set_xlabel('Offset [AU]', fontsize=14)
# ax[1].set_xticks([0, conv_image.shape[0]//2, conv_image.shape[0]-1])
# ax[1].set_xticklabels([-np.round(sizeau//2), 0, np.round(sizeau//2)])
# # ax[1].set_ylabel('Offset [AU]', fontsize=14)
# # ax[1].text(0.9, 0.9, 'eDisk (1.3 mm)', transform=ax[1].transAxes, fontsize=14, color='black')
# ax[1].set_title('Model (Accretion)', fontsize=14)


# residual_ax = ax[2].imshow(residual*1e3, origin='lower', cmap = residual_cmp, vmin=-2, vmax=2)
# colorbar = fig.colorbar(residual_ax, ax=ax[2], pad=0.00, aspect=30, shrink=.98)
# colorbar.set_label('Intensity (mJy/beam)')
# ax[2].set_xlabel('Offset [AU]', fontsize=14)
# ax[2].set_xticks([0, residual.shape[0]//2, residual.shape[0]-1])
# ax[2].set_xticklabels([-np.round(sizeau//2), 0, np.round(sizeau//2)])
# # ax[2].set_ylabel('Offset [AU]', fontsize=14)
# # ax[2].text(0.7, 0.7, 'Residual', transform=ax[2].transAxes, fontsize=14, color='black')
# ax[2].set_title('Residual', fontsize=14)
# plt.savefig('residual_accretion_1e-5.pdf',transparent=True)
# plt.close("all")








# a = 0.1
# L_star = 3.0
# Q = 0.3
# mdot = 1e-6
# heat = 'radiation'


# # model = generate_model(a, L_star, Q, mdot, heat = heat)
# model_im = image.readImage(fname=f'./simulation/outfile/conti_a_{a}_Lstar_{L_star}_Q_{Q}_mdot_{mdot}_{heat}_scat.out')
# npix = model_im.nx
# pixel_area = (sizeau/npix/140)**2
# conv_image = model_im.imConv(dpc=140, fwhm=beam_axis, pa=-beam_pa)
# conv_image = conv_image.imageJyppix * beam_area/pixel_area/(140**2)

# f_interp_model = interp2d(np.linspace(-(crop_sizeau//2), crop_sizeau//2, conv_image.shape[0]),
#                             np.linspace(-(crop_sizeau//2), crop_sizeau//2, conv_image.shape[1]), conv_image, kind='linear')
# conv_image = f_interp_model(np.linspace(-(crop_sizeau//2), crop_sizeau//2, 1000),
#                             np.linspace(-(crop_sizeau//2), crop_sizeau//2, 1000))
# peak_idx_x_model, peak_idx_y_model = np.unravel_index(np.argmax(conv_image, axis=None), conv_image.shape)
# conv_image = conv_image[peak_idx_x_model-400:peak_idx_x_model+400, peak_idx_y_model-400:peak_idx_y_model+400].T
# mask = conv_image < 5*sigma
# conv_image[mask] = np.nan

# residual = edisk_image - conv_image





# fig, ax = plt.subplots(1,3, sharex=False, sharey=True, figsize=(15,5))
# fig.subplots_adjust(left=0.05, right=0.97, top=0.9, bottom=0.1, wspace=0.0, hspace=0.0)


# edisk = ax[0].imshow(edisk_image*1e3, origin='lower', cmap='plasma', vmin=0, vmax=5)
# colorbar = fig.colorbar(edisk, ax=ax[0], pad=0.00, aspect=30, shrink=.98)
# # colorbar.set_label('Intensity (mJy/beam)')
# ax[0].set_xticks([0, edisk_image.shape[0]//2, edisk_image.shape[0]-1])
# ax[0].set_xticklabels([-np.round(sizeau//2), 0, np.round(sizeau//2)])
# ax[0].set_xlabel('Offset [AU]', fontsize=14)
# ax[0].set_yticks([0, edisk_image.shape[0]//2, edisk_image.shape[0]-1])
# ax[0].set_yticklabels([-np.round(sizeau//2), 0, np.round(sizeau//2)])
# ax[0].set_ylabel('Offset [AU]', fontsize=14)
# # ax[0].text(0.9, 0.9, 'eDisk (1.3 mm)', transform=ax[0].transAxes, fontsize=14, color='black')
# ax[0].set_title('CB68', fontsize=14)

# model = ax[1].imshow(conv_image*1e3, origin='lower', cmap='plasma', vmin=0, vmax=5)
# colorbar = fig.colorbar(model, ax=ax[1], pad=0.00, aspect=30, shrink=.98)
# # colorbar.set_label('Intensity (mJy/beam)')
# ax[1].set_xlabel('Offset [AU]', fontsize=14)
# ax[1].set_xticks([0, conv_image.shape[0]//2, conv_image.shape[0]-1])
# ax[1].set_xticklabels([-np.round(sizeau//2), 0, np.round(sizeau//2)])
# # ax[1].set_ylabel('Offset [AU]', fontsize=14)
# # ax[1].text(0.9, 0.9, 'eDisk (1.3 mm)', transform=ax[1].transAxes, fontsize=14, color='black')
# ax[1].set_title('Model (Radiation)', fontsize=14)


# residual_ax = ax[2].imshow(residual*1e3, origin='lower', cmap = residual_cmp, vmin=-2, vmax=2)
# colorbar = fig.colorbar(residual_ax, ax=ax[2], pad=0.00, aspect=30, shrink=.98)
# colorbar.set_label('Intensity (mJy/beam)')
# ax[2].set_xlabel('Offset [AU]', fontsize=14)
# ax[2].set_xticks([0, residual.shape[0]//2, residual.shape[0]-1])
# ax[2].set_xticklabels([-np.round(sizeau//2), 0, np.round(sizeau//2)])
# # ax[2].set_ylabel('Offset [AU]', fontsize=14)
# # ax[2].text(0.7, 0.7, 'Residual', transform=ax[2].transAxes, fontsize=14, color='black')
# ax[2].set_title('Residual', fontsize=14)
# plt.savefig('residual_radiation_Q_0.3.pdf',transparent=True)
# plt.close("all")