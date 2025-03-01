import numpy as np
import matplotlib.pyplot as plt
from generate_cube import *

from CB68.data_dict import data_dict

# methanol_data = data_dict["ch3oh_218_faust"]
# beam_axis = [methanol_data["bmaj"], methanol_data["bmin"]]
# beam_pa = methanol_data["bpa"]

# sizeau = 280
# npix   = 500
# pixel_area = (sizeau/npix/140)**2
# beam_area = beam_axis[0]*beam_axis[1]*np.pi/(4*np.log(2))

# a_list = [1e1, 5e0, 1e0, 5e-1, 1e-1, 5e-2, 1e-2]
# mdot_list = [1e-6, 5e-7, 1e-7, 5e-8, 1e-8]
# Q_list = [2, 1.5, 1, 0.5]
# snowline_list = [100, None]

# for a in a_list:
#     for mdot in mdot_list:
#         for Q in Q_list:
#             for snowline in snowline_list:
#                 generate_model( amax     =  a, # mm
#                                 mstar    =  0.14, # Msun
#                                 mdot     =  mdot, # Msun/yr
#                                 Q        =  Q, # Toomre Q
#                                 snowline =  snowline, # temperature of sublimation
#                                 rcb      =  None)

#                 produce_cube(fname=f"a_{a}_mdot_{mdot}_Q_{Q}_snowline_{snowline}",channel=False, pv=True)


# methanol_data = data_dict["ch3oh_218_faust"]
beam_axis = [0.497, 0.385]
beam_pa = -62.5

def convolve(fname):
  with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
    im = image.readImage(fname=fname)
  im_conv = im.imConv(dpc=140, fwhm=beam_axis, pa=-beam_pa)
  pixel_area = np.abs(model_image.sizepix_x*model_image.sizepix_y)
  return im, im_conv

def rotate_image(image, posang):
  if isinstance(image, np.ndarray)!=True:
    image.imageJyppix= ndimage.rotate(image.imageJyppix, posang, reshape=False, axes=(1, 0))
    image.imageJyppix = np.nan_to_num(image.imageJyppix, nan=0)
    return image.imageJyppix
  else:
    image = ndimage.rotate(image, posang, reshape=False, axes=(1, 0))
    image = np.nan_to_num(image, nan=0)
    return image

def pv(image_array, average_width=3):
  center = image_array.shape[1]//2
  pv_slice = np.mean(image_array[:, center-average_width//2:center+average_width//2, :], axis=1)
  return pv_slice

def write_fits(image, fname):
  if os.path.exists(fname):
    os.remove(fname)
  image.writeFits(fname=fname, dpc=140, coord='16h57m19.647s -16d09m23.94s')

def initialize_image(fname, convolve=True):
  if convolve:
    im, im_conv = convolve(fname)
  else:
    with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
      im = image.readImage(fname=fname)
    im_conv = im
  return im, im_conv



sizeau = 280
npix   = 500
pixel_area = (sizeau/npix/140)**2
beam_area = beam_axis[0]*beam_axis[1]*np.pi/(4*np.log(2))

a_list = [1e1, 5e0, 1e0, 5e-1, 1e-1, 5e-2, 1e-2]
mdot_list = [1e-6, 5e-7, 1e-7, 5e-8, 1e-8]
Q_list = [2, 1.5, 1, 0.5]
snowline_list = [100, None]

for a in a_list:
    for mdot in mdot_list:
        for Q in Q_list:
            for snowline in snowline_list:
              try:
                im_line = image.readImage(fname=f'./test/outfile/pv_a_{a}_mdot_{mdot}_Q_{Q}_snowline_{snowline}_scat.out')
                model_image_line = im_line.imConv(dpc=140, fwhm=beam_axis, pa=-beam_pa)

                im_conti = image.readImage(fname=f'./test/outfile/pv_a_{a}_mdot_{mdot}_Q_{Q}_snowline_{snowline}_conti.out')
                model_image_conti = im_conti.imConv(dpc=140, fwhm=beam_axis, pa=-beam_pa)
                model_image_conti.imageJyppix = np.tile(model_image_conti.imageJyppix, (1, 1, im_line.nwav))

                model_image = model_image_line.imageJyppix - model_image_conti.imageJyppix
                model_image = model_image/(140**2)
                rotated_image = rotate_image(model_image, 45)
                x_axis = np.linspace(-2, 2, rotated_image.shape[0])
                v_axis = np.linspace(-5, 15, im_line.nwav)
                pv_slice = pv(rotated_image, 5)
                plt.pcolormesh(x_axis, v_axis, pv_slice.T, 
                              shading="nearest", rasterized=True, cmap='gist_ncar',
                              vmin=0, vmax=0.04)
                plt.xlabel("Offset [au]",fontsize = 16)
                plt.ylabel("Velocity [km/s]",fontsize = 16)
                plt.plot([0, 0], plt.ylim(), 'w:')
                plt.plot(plt.xlim(), [5, 5], 'w:')
                plt.colorbar(label='Jy/pixel')
                plt.savefig(f'./figures/a_{a}_mdot_{mdot}_Q_{Q}_snowline_{snowline}.pdf', transparent=True)
                plt.close("all")
              except:
                pass


