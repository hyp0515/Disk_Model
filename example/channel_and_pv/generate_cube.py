import numpy as np
import matplotlib.pyplot as plt
from scipy import ndimage
import sys
import os
import io
import contextlib
sys.path.insert(0,'../../')
from radmc.setup import radmc3d_setup
from radmc.simulate import generate_simulation
from radmc3dPy import image
pc  = 3.08572e18
from CB68.data_dict import data_dict
import astropy.constants as const
cc = const.c.cgs.value

class general_parameters:
    '''
    A class to store the parameters for individual kinds of grids.
    Details of individual parameters should refer to the functions that generate the grids.
    '''
    def __init__(self, **kwargs
                 ):
        for k, v in kwargs.items():
          # add parameters as attributes of this object
          setattr(self, k, v)

    def __del__(self):
      pass

    def add_attributes(self, **kwargs):
      '''
      Use this function to set the values of the attributs n1, n2, n3,
      which are number of pixels in the first, second, and third axes. 
      '''
      for k, v in kwargs.items():
        # add parameters as attributes of this object
        setattr(self, k, v)


def generate_model(amax=0.05, l_star=5.0, Q=0.5, mdot = 1e-7, heat="radiation"):

  
  model = radmc3d_setup(silent=False)
  model.get_mastercontrol(filename=None,
                          comment=None,
                          incl_dust=1,
                          incl_lines=1,
                          nphot=1000000,
                          nphot_scat=10000000,
                          scattering_mode_max=2,
                          istar_sphere=1,
                          num_cpu=None,
                          modified_random_walk = 1
                          )
  model.get_linecontrol(filename=None,
                      methanol='ch3oh leiden 0 0 0')
  model.get_continuumlambda(filename=None,
                          comment=None,
                          lambda_micron=None,
                          append=False)

  model.get_diskcontrol(  d_to_g_ratio = 0.01,
                          a_max=amax, 
                          Mass_of_star=0.14, 
                          Accretion_rate=mdot,
                          Radius_of_disk=25,
                          NR=200,
                          NTheta=200,
                          NPhi=20,
                          Q=Q)
  model.get_vfieldcontrol(Kep=True,
                          vinfall=0.5,
                          Rcb=None,
                          outflow=None)
  model.get_heatcontrol(L_star=l_star,
                        R_star=1,
                        heat=heat)
  model.get_gasdensitycontrol(abundance=1e-10,
                              snowline=100,
                              enhancement=1e5,
                              gas_inside_rcb=True)

##############################################

def produce_cube(fname, channel=True, pv=False):
  condition_parms = general_parameters(
      nodust      = False,
      scat        = True,
      extract_gas = True,
  )


  simulate_mutual_parms = {
      "incl"      : 73,
      "line"      : 240,
      "npix"      : 500,
      "sizeau"    : 280,  # 2 arcsec with 140 pc distance
      "v_width"   : 5,
      "vkms"      : 0,
      "posang"    : 45,
      "dir"       : './simulation/',
      "fname"     : fname,
  }

  channel_cube_parms = general_parameters(
      **simulate_mutual_parms,
      nlam=11,
  )

  pv_cube_parms = general_parameters(
      **simulate_mutual_parms,
      nlam=50,
  )

  simulation_parms = general_parameters(
      condition_parms    = condition_parms,
      channel_cube_parms = channel_cube_parms,
      pv_cube_parms      = pv_cube_parms,
      save_out=True,
      save_npz=True,
  )

  simulation = generate_simulation(
      parms=simulation_parms,
      channel       = channel,
      pv            = pv,
  )

##############################################

# generate_model(amax=0.05, l_star=5.0, Q=0.5, mdot = 1e-7, heat="radiation")
# produce_cube(fname="irr_best", channel=True, pv=True)

# generate_model(amax=10.0, l_star=0.1, Q=0.5, mdot = 1e-6, heat="accretion")
# produce_cube(fname="acc_best", channel=True, pv=True)

# a_list = [5e-1, 1e0, 5e0, 1e1]
# for a in a_list:
#   generate_model(amax=a, l_star=5.0, Q=0.5, mdot = 1e-7, heat="radiation")
#   produce_cube(fname=f"irr_{a}", channel=True, pv=True)
# a_list = [1e-2, 5e-1, 1e-1, 5e-1, 1e0, 5e0] 
# for a in a_list:
#   generate_model(amax=a, l_star=0.1, Q=0.5, mdot = 1e-6, heat="accretion")
#   produce_cube(fname=f"acc_{a}", channel=True, pv=True)
# ##############################################
# methanol_data = data_dict["ch3oh_218_faust"]
# beam_axis = [methanol_data["bmaj"], methanol_data["bmin"]]
# beam_pa = methanol_data["bpa"]

# def convolve(fname):
#   with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
#     im = image.readImage(fname=fname)
#   im_conv = im.imConv(dpc=140, fwhm=beam_axis, pa=-beam_pa)
#   pixel_area = np.abs(model_image.sizepix_x*model_image.sizepix_y)
#   return im, im_conv

# def rotate_image(image, posang):
#   if isinstance(image, np.ndarray)!=True:
#     image.imageJyppix= ndimage.rotate(image.imageJyppix, posang, reshape=False, axes=(1, 0))
#     image.imageJyppix = np.nan_to_num(image.imageJyppix, nan=0)
#     return image.imageJyppix
#   else:
#     image = ndimage.rotate(image, posang, reshape=False, axes=(1, 0))
#     image = np.nan_to_num(image, nan=0)
#     return image

# def pv(image_array, average_width=10):
#   center = image_array.shape[1]//2
#   pv_slice = np.mean(image_array[:, center-average_width//2:center+average_width//2, :], axis=1)
#   return pv_slice

# # def write_fits(image, fname):
# #   if os.path.exists(fname):
# #     os.remove(fname)
# #   image.writeFits(fname=fname, dpc=140, coord='16h57m19.647s -16d09m23.94s')

# def initialize_image(fname, convolve=True):
#   if convolve:
#     im, im_conv = convolve(fname)
#   else:
#     with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
#       im = image.readImage(fname=fname)
#     im_conv = im
#   return im, im_conv



# sizeau = 280
# npix   = 500
# pixel_area = (sizeau/npix/140)**2
# beam_area = beam_axis[0]*beam_axis[1]*np.pi/(4*np.log(2))



# im_line = image.readImage(fname='./simulation/outfile/channel_irr_best_scat.out')
# model_image_line = im_line.imConv(dpc=140, fwhm=beam_axis, pa=-beam_pa)

# im_conti = image.readImage(fname='./simulation/outfile/channel_irr_best_conti.out')
# # im_conti.imageJyppix = np.tile(im_conti.imageJyppix, (1, 1, im_line.nwav))
# model_image_conti = im_conti.imConv(dpc=140, fwhm=beam_axis, pa=-beam_pa)
# model_image_conti.imageJyppix = np.tile(model_image_conti.imageJyppix, (1, 1, im_line.nwav))

# model_image = model_image_line.imageJyppix - model_image_conti.imageJyppix
# model_image = model_image/(140**2)*beam_area/pixel_area
# rotated_image = rotate_image(model_image, 45)
# # plt.imshow(rotated_image[:, :, 5].T, origin='lower')
# # plt.show()

# x_axis = np.linspace(-1, 1, rotated_image.shape[0])
# v_axis = np.linspace(-5, 15, im_line.nwav)
# pv_slice = pv(rotated_image, 5)
# plt.pcolormesh(x_axis, v_axis, pv_slice.T, 
#                shading="nearest", rasterized=True, cmap='gist_ncar',
#                vmin=0, vmax=0.02)
# plt.xlabel("Offset [au]",fontsize = 16)
# plt.ylabel("Velocity [km/s]",fontsize = 16)
# plt.plot([0, 0], plt.ylim(), 'w:')
# plt.plot(plt.xlim(), [5, 5], 'w:')
# plt.show()
# print(im_line.freq)
# print(cc / 1e5 * (2.18440063e11 - im_line.freq) / 2.18440063e11)
# fig, ax = plt.subplots(1, im_line.nwav)
# for i in range(im_line.nwav):
#   ax[i].imshow(model_image[:, :, i].T, origin='lower', vmin=0, vmax=0.01)
# plt.show()

# fig, ax = plt.subplots(1, 2)
# ax[0].imshow(model_image.imageJyppix[:, :, 0].T, origin='lower')
# rotated_image = rotate_image(model_image, 45)
# ax[1].imshow(rotated_image[:, :, 0].T, origin='lower')
# plt.show()
# # pv_slice = pv(rotated_image, 3)


# plt.pcolormesh(pv_slice, origin='lower')
# plt.show()
# model_image.image = model_image.image/(140**2)
# im = convolve('./test/outfile/channel_test_scat.out')
# write_fits(im, './channel_test_scat_noconv.fits')