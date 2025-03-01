import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.insert(0,'../../')
from radmc.setup import radmc3d_setup
from radmc.simulate import generate_simulation
from radmc.plot import generate_plot
from radmc3dPy import image
import os 
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


model = radmc3d_setup(silent=False)
model.get_mastercontrol(filename=None,
                        comment=None,
                        incl_dust=1,
                        incl_lines=1,
                        nphot=100000,
                        nphot_scat=5000000,
                        scattering_mode_max=2,
                        istar_sphere=1,
                        num_cpu=None,
                        modified_random_walk = 1)
model.get_linecontrol(filename=None,
                      methanol='ch3oh leiden 0 0 0')
model.get_continuumlambda(filename=None,
                          comment=None,
                          lambda_micron=None,
                          append=False)

model.get_diskcontrol(  d_to_g_ratio = 0.01,
                        a_max=.01, 
                        Mass_of_star=0.14, 
                        Accretion_rate=1e-6,
                        Radius_of_disk=30,
                        NR=200,
                        NTheta=200,
                        NPhi=10)
model.get_vfieldcontrol(Kep=True,
                        vinfall=0.5,
                        Rcb=None,
                        outflow=None)
model.get_heatcontrol(heat='accretion')
model.get_gasdensitycontrol(abundance=1e-10,
                            snowline=100,
                            enhancement=1e5,
                            gas_inside_rcb=True)

# os.system('radmc3d mctherm')

##############################################

condition_parms = general_parameters(
    nodust      = False,
    scat        = True,
    extract_gas = True,
)


simulate_mutual_parms = {
    "incl"      : 70,
    "line"      : 240,
    "npix"      : 500,
    "sizeau"    : 1000,
    "v_width"   : 10,
    "vkms"      : 0,
    "v_width"   : 10,
    "dir"       : './test/',
    "fname"     : 'test',
}

channel_cube_parms = general_parameters(
    **simulate_mutual_parms,
    nlam=11,
)

pv_cube_parms = general_parameters(
    **simulate_mutual_parms,
    nlam=50,
)

sed_parms = general_parameters(
    **simulate_mutual_parms, 
    scat=True,
    freq_min=5e1, freq_max=5e2, nlam=10,
)

spectrum_parms = general_parameters(
    **simulate_mutual_parms,
    nlam=10
)
conti_parms = general_parameters(
    **simulate_mutual_parms,
    wav=1300,
    scat=True,
)

simulation_parms = general_parameters(
    condition_parms    = condition_parms,
    channel_cube_parms = channel_cube_parms,
    pv_cube_parms      = pv_cube_parms,
    conti_parms        = conti_parms,
    sed_parms          = sed_parms,
    spectrum_parms     = spectrum_parms,
    save_out=True,
    save_npz=True,
)

simulation = generate_simulation(
    parms=simulation_parms,
    channel       = False,
    pv            = False,
    conti         = True,
    sed           = False,
    line_spectrum = False
)
##############################################


##############################################
# plot_mutual_parms = {
   
# }


# plot_profile_parms = general_parameters(
#     dir=f'./test/', fname='profile'
# )

# plot_channel_parms = general_parameters(
#     cube_dir = './test/npzfile/',
#     cube_fname = [f'channel_test_scat.npz',
#                     f'channel_test_conti.npz'],
#     extracted = True,
#     conti = True,
#     convolve = True,
#     vkms = 5,
#     fwhm = 60,
#     d = 140,
#     title = f'test',
#     dir = f'./test/test/',
#     fname = 'channel'
# )

# plot_pv_parms = general_parameters(
#     cube_dir = './test/npzfile/',
#     cube_fname = [f'pv_test_scat.npz',
#                     f'pv_test_conti.npz'],
#     extracted = True,
#     conti = True,
#     convolve = True,
#     vkms = 5,
#     fwhm = 60,
#     d = 140,
#     CB68 = True,
#     title = f'test',
#     dir = f'./test/',
#     fname = 'pv'
# )

# plot_sed_parms = general_parameters(
#     sed_dir = './test/npzfile/',
#     sed_fname = f'sed_test_scat.npz',
#     d = 140,
#     CB68 = True,
#     title = f'test',
#     dir = f'./test/',
#     fname = 'sed'
# )

# plot_spectrum_parms = general_parameters(
#     spectra_dir = './test/npzfile/',
#     spectra_fname = [f'spectrum_test_scat.npz',
#                     f'spectrum_test_conti.npz'],
#     extracted = True,
#     vkms = 5,
#     d = 140,
#     title = f'test',
#     dir = f'./test/',
#     fname = 'spectrum'
# )

# plot_continuum_parms = general_parameters(
#     conti_dir = './test/npzfile/',
#     conti_fname = f'conti_test_scat.npz',
#     convolve = True,
#     fwhm = 60,
#     d = 140,
#     title = f'test',
#     dir = f'./test/',
#     fname = 'conti'
# )


# plot_parms = general_parameters(
#     profile_parms   = plot_profile_parms,
#     channel_parms   = plot_channel_parms,
#     pv_parms        = plot_pv_parms,
#     sed_parms       = plot_sed_parms,
#     spectra_parms   = plot_spectrum_parms,
#     continuum_parms = plot_continuum_parms
# )

# plot = generate_plot(
#     parms     = plot_parms, 
#     profile   = False,
#     channel   = True,
#     pv        = False,
#     continuum = False,
#     sed       = False,
#     spectrum  = False,
# )

