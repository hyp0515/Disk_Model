import numpy as np
import matplotlib.pyplot as plt
from scipy import ndimage
from radmc3dPy import image
# from CB68.data_dict import data_dict
from make_conti import generate_model

def rotate_image(image, posang):
  if isinstance(image, np.ndarray)!=True:
    image.imageJyppix= ndimage.rotate(image.imageJyppix, posang, reshape=False, axes=(1, 0))
    image.imageJyppix = np.nan_to_num(image.imageJyppix, nan=0)
    return image.imageJyppix
  else:
    image = ndimage.rotate(image, posang, reshape=False, axes=(1, 0))
    image = np.nan_to_num(image, nan=0)
    return image

def radial_intensity(image_array, center, width):
  radial_profile = np.mean(image_array[:, center-width//2:center+width//2], axis=1)
  return radial_profile

sizeau = 80
npix = 500
pixel_area = (sizeau/npix/140)**2
beam_axis = [0.0363, 0.0274]
beam_area = beam_axis[0]*beam_axis[1]*np.pi/(4*np.log(2))
r_axis = np.linspace(-(sizeau//2), sizeau//2, npix, endpoint=True)

edisk_radial = np.load("edisk_radial.npz")


heat_list = ["accretion", "radiation", "combine"]
a_list    = [5e1 ,1e1, 5e0, 1e0, 5e-1, 1e-1, 5e-2, 1e-2]
for a in a_list:
    plt.figure(figsize=(8, 4))
    for heat in heat_list:
        generate_model(amax=a, heat=heat)
        im = image.readImage(fname=f'./test/outfile/conti_a_{a}_{heat}_scat.out')
        im_conv = im.imConv(dpc=140, fwhm=beam_axis, pa=-69.4)
        im_conv = rotate_image(im_conv, 45)
        i_r = radial_intensity(im_conv, npix//2, 10)
        plt.plot(r_axis, i_r*1e3*(beam_area/pixel_area)/(140**2), label=f"{heat}")
    plt.plot(edisk_radial["r_axis"], edisk_radial["i_r"]*1e3, label="eDisk", linestyle="--")
    plt.xlim((-(sizeau//2), sizeau//2))
    plt.ylim(bottom=0)
    plt.xlabel("Offset (AU)")
    plt.ylabel("Intensity (mJy/beam)")
    plt.title(f"Continuum radial profile along the major axis "+r"($a_{max}=$"+f"{a}mm)")
    plt.legend()
    plt.savefig(f'./figures/a_{a}.pdf', transparent=True)
    plt.close("all")


# generate_model(amax=0.3, heat="irradiation")
# im = image.readImage(fname=f'./test/outfile/conti_a_0.3_irradiation_scat.out')
# im_conv = im.imConv(dpc=140, fwhm=beam_axis, pa=-69.4)
# im_conv = rotate_image(im_conv, 45)
# # im_conv[:,npix//2-5//2:npix//2+5//2,0] = np.tile(np.linspace(0, 0.01, npix)[:,np.newaxis], (1, 4))
# plt.imshow(im_conv[:,:,0].T, origin="lower", cmap="inferno")
# plt.show()