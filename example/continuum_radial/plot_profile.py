import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
import matplotlib.font_manager as fm
from scipy import ndimage
from scipy.optimize import curve_fit
import os
import sys
from scipy.interpolate import interp1d
import corner
from radmc3dPy import image
from radmc3dPy.analyze import *
# from CB68.data_dict import data_dict
from make_conti import generate_model
sys.path.append("..")
from profile_radial.make_disk import generate_disk




sizeau = 80
npix = 500
pixel_area = (sizeau/npix/140)**2
beam_axis = [0.0363, 0.0274]
beam_area = beam_axis[0]*beam_axis[1]*np.pi/(4*np.log(2))
r_axis = np.linspace(-(sizeau//2), sizeau//2, npix, endpoint=True)


edisk_radial = np.load("edisk_radial.npz")
i_r_obs = edisk_radial["i_r"]



def rotate_image(image, posang):
  if isinstance(image, np.ndarray)!=True:
    image.imageJyppix= ndimage.rotate(image.imageJyppix, posang, reshape=False, axes=(1, 0))
    image.imageJyppix = np.nan_to_num(image.imageJyppix, nan=0)
    return image.imageJyppix[:,:,0]
  else:
    image = ndimage.rotate(image, posang, reshape=False, axes=(1, 0))
    image = np.nan_to_num(image, nan=0)
    return image

def radial_intensity(image_array, center, width):
  if center is None:
    peak_idx_x, peak_idx_y = np.unravel_index(np.argmax(image_array, axis=None), image_array.shape)
    center = peak_idx_y
  radial_profile = np.mean(image_array[:, center-width//2:center+width//2], axis=1)
  return radial_profile

def get_radial_profile(image, beam_axis, posang=45, center=None, width=10):
  conv_image = image.imConv(dpc=140, fwhm=beam_axis, pa=-69.4)
  conv_image = rotate_image(conv_image, posang)
  conv_image *= beam_area/pixel_area/(140**2)
  i_r = radial_intensity(conv_image, center, width)
  return i_r

def chi(i_r_model, i_r_obs):
    if len(i_r_model) != len(i_r_obs):
        interp_func = interp1d(np.linspace(0, 1, len(i_r_model)), i_r_model, kind='cubic')
        i_r_model = interp_func(np.linspace(0, 1, len(i_r_obs)))
    return np.sum(((i_r_model - i_r_obs)**2)/(21e-6**2))

def plot_i_r(i_r_model, a, L_star, Q, mdot, heat='radiation'):
    plt.plot(r_axis, i_r_model*1e3, label=" radiation model")
    plt.plot(edisk_radial["r_axis"], edisk_radial["i_r"]*1e3, label="eDisk", linestyle="--")
    plt.xlim((-(sizeau//2), sizeau//2))
    plt.ylim(bottom=0)
    plt.xlabel("Offset (AU)")
    plt.ylabel("Intensity (mJy/beam)")
    plt.title(f"Continuum radial profile along the major axis")
    plt.legend()
    plt.savefig(f"./figures/i_r/{heat}/a_{a}_Lstar_{L_star}_Q_{Q}_mdot_{mdot}_{heat}.pdf", transparent=True)
    plt.close("all")

def plot_image(image, beam_axis, posang, a, L_star, Q, mdot, heat='radiation'):
    conv_image = image.imConv(dpc=140, fwhm=beam_axis, pa=-69.4)
    conv_image = rotate_image(conv_image, posang)
    conv_image *= beam_area/pixel_area/(140**2)
    plt.imshow(conv_image.T, origin='lower', cmap="inferno", vmin=10*21e-6)
    plt.colorbar()
    plt.savefig(f"./figures/image/{heat}/a_{a}_Lstar_{L_star}_Q_{Q}_mdot_{mdot}_{heat}.pdf", transparent=True)
    plt.close("all")


def read_from_dir(dir):
    os.chdir(dir)
    d = readData(dtemp=True, ddens=True, gdens=True, ispec='ch3oh')
    grid = readGrid(wgrid=False)
    os.chdir('..')
    return d, grid

def plot_profile(d, grid, d2=None):
    nch3oh    = d.ndens_mol[:, :, 0, 0]
    dust      = np.sum(d.rhodust[:, :, 0, :], axis=2)
    t         = np.mean(d.dusttemp[:, :, 0, :], axis=2)
    if d2 is not None:
        nch3oh2    = d2.ndens_mol[:, :, 0, 0]
        dust2      = np.sum(d2.rhodust[:, :, 0, :], axis=2)
        t2         = np.mean(d2.dusttemp[:, :, 0, :], axis=2)
    R, Theta, Phi =  grid.x/au, grid.y, grid.z
    fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(18, 6),
                        subplot_kw={'projection': 'polar'})
    fig.subplots_adjust(left=0.05, right=0.95, top=0.9, bottom=0.1, wspace=0.1, hspace=0.05)
    cmaps = ['BuPu', 'OrRd', 'BuPu']
    titles = [r'$\rho_{dust}$', r'$T$', r'$n_{\mathregular{CH_3OH}}$']
    cbar = [r'log($\rho$) [g$cm^{-3}$]', r'log(T) [K]',r'log($n_{\mathregular{CH_3OH}}$) [$cm^{-3}$]']

    if d2 is None:
        for idx_val, val in enumerate([dust, t, nch3oh]):
            c = ax[idx_val].pcolormesh(Theta-np.pi/2, R, np.log10(val), shading='auto', cmap=cmaps[idx_val])
            ax[idx_val].pcolormesh(Theta+np.pi/2, R, np.log10(val), shading='auto', cmap=cmaps[idx_val])
            if idx_val == 0:
                den = val
                # levels = np.linspace(np.log10(den).min(), np.log10(den).max(), 3)
                levels = [-20, -15, -10]
            # ax[idx_val].contour(Theta-np.pi/2, R, np.log10(den), levels=levels, colors='k', linewidths=.7, linestyles='dashed')
            # ax[idx_val].contour(Theta+np.pi/2, R, np.log10(den), levels=levels, colors='k', linewidths=.7, linestyles='dashed')
            ax[idx_val].set_xticks([])
            ax[idx_val].set_yticks([])
            fig.colorbar(c, ax=ax[idx_val], orientation='vertical', shrink=0.7).set_label(cbar[idx_val], fontsize=18)
            ax[idx_val].set_title(titles[idx_val], fontsize=26, color='k')
    else:
        for idx_val, val in enumerate([(dust, dust2), (t, t2), (nch3oh, nch3oh2)]):
            val1, val2 = val
            # if idx_val == 1:
            #     c = ax[idx_val].pcolormesh(Theta-np.pi/2, R, val1, shading='auto', cmap=cmaps[idx_val])
            #     ax[idx_val].pcolormesh(Theta+np.pi/2, R, val2, shading='auto', cmap=cmaps[idx_val])
            # else:
            c = ax[idx_val].pcolormesh(Theta-np.pi/2, R, np.log10(val1), shading='auto', cmap=cmaps[idx_val])
            ax[idx_val].pcolormesh(Theta+np.pi/2, R, np.log10(val2), shading='auto', cmap=cmaps[idx_val])
            if idx_val == 0:
                den = val1
                # levels = np.linspace(np.log10(den).min(), np.log10(den).max(), 3)
                levels = [-20, -15, -12]
            # ax[idx_val].contour(Theta-np.pi/2, R, np.log10(den), levels=levels, colors='k', linewidths=.5, linestyles='dotted')
            # ax[idx_val].contour(Theta+np.pi/2, R, np.log10(den), levels=levels, colors='k', linewidths=.5, linestyles='dotted')
            ax[idx_val].set_xticks([])
            ax[idx_val].set_yticks([])
            fig.colorbar(c, ax=ax[idx_val], orientation='vertical', shrink=0.7).set_label(cbar[idx_val], fontsize=18)
            ax[idx_val].set_title(titles[idx_val], fontsize=26, color='k')


    scale_bar_ax = fig.add_axes([0.48, 0.10, 0.1, 0.02]) # [left, bottom, width, height]
    scale_bar = AnchoredSizeBar(scale_bar_ax.transData,
                                1,  # Size of the scale bar in data coordinates
                                f'{round(R[-1])} AU',  # Label for the scale bar
                                'lower center',  # Location
                                pad=0.1,
                                color='black',
                                frameon=False,
                                size_vertical=0.01,
                                fontproperties=fm.FontProperties(size=12))

    scale_bar_ax.add_artist(scale_bar)
    scale_bar_ax.set_axis_off()

    scale_bar_ax = fig.add_axes([0.17, 0.10, 0.1, 0.02]) # [left, bottom, width, height]
    scale_bar = AnchoredSizeBar(scale_bar_ax.transData,
                                1,  # Size of the scale bar in data coordinates
                                f'{round(R[-1])} AU',  # Label for the scale bar
                                'lower center',  # Location
                                pad=0.1,
                                color='black',
                                frameon=False,
                                size_vertical=0.01,
                                fontproperties=fm.FontProperties(size=12))

    scale_bar_ax.add_artist(scale_bar)
    scale_bar_ax.set_axis_off()

    scale_bar_ax = fig.add_axes([0.78, 0.10, 0.1, 0.02]) # [left, bottom, width, height]
    scale_bar = AnchoredSizeBar(scale_bar_ax.transData,
                                1,  # Size of the scale bar in data coordinates
                                f'{round(R[-1])} AU',  # Label for the scale bar
                                'lower center',  # Location
                                pad=0.1,
                                color='black',
                                frameon=False,
                                size_vertical=0.01,
                                fontproperties=fm.FontProperties(size=12))

    scale_bar_ax.add_artist(scale_bar)
    scale_bar_ax.set_axis_off()

def gaussian(r, I0, r0, sigma):
    return I0 * np.exp(-((r - r0) ** 2) / (2 * sigma ** 2))

def gaussian_fit(i_r, r_axis, extract_index=None):
    I0_guess = np.max(i_r)
    r0_guess = r_axis[np.argmax(i_r)]
    sigma_guess = (np.max(r_axis) - np.min(r_axis)) / 4  # Rough estimate
    p0 = [I0_guess, r0_guess, sigma_guess]
    if extract_index is not None:
        i_r = i_r[extract_index:-extract_index]
        r_axis = r_axis[extract_index:-extract_index]
    popt, pcov = curve_fit(gaussian, r_axis, i_r, p0=p0)
    return popt


a = 0.1
L_star = 3
Q = 1.5
mdot = 1e-6
heat = "radiation"


# disk_acc = generate_disk(
#     amax   = 0.05,
#     mstar  = 0.14,
#     mdot   = 1e-7,
#     rd     = 25,
#     Q      = 0.5,
#     l_star = 5.0,
#     r_star = 1,
#     heat   = "accretion",
#     dir="./disk_acc/"
# )

# disk_irr = generate_disk(
#     amax   = 0.05,
#     mstar  = 0.14,
#     mdot   = 1e-7,
#     rd     = 25,
#     Q      = 0.5,
#     l_star = 5.0,
#     r_star = 1,
#     heat   = "irradiation",
#     dir="./disk_irr/"
# )

d_acc, grid = read_from_dir("./disk_acc/")
d_irr, _    = read_from_dir("./disk_irr/")
plot_profile(d_acc, grid, d_irr)

# plt.grid(False)
plt.savefig("disk_profile.pdf", transparent=True, bbox_inches='tight', dpi=100)


