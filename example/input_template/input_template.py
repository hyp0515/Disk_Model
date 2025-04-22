import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.insert(0,'../../')
from radmc.setup import radmc3d_setup
from radmc.simulate import generate_simulation
from radmc.plot import generate_plot
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
                        a_max           = 0.01, # mm
                        Mass_of_star    = 0.14, # Msun
                        Accretion_rate  = 1e-7, # Msun/yr
                        Radius_of_disk  = 25,   # AU
                        NR    =200,
                        NTheta=200,
                        NPhi  =20,
                        Q=1.5 # Toomre Q
                        )
model.get_vfieldcontrol(Kep=True,
                        vinfall=0.5, # the infall velocity (unit: Keperian velocity)
                        Rcb=None, # the centrifugal barrier
                        outflow=None)
model.get_heatcontrol(L_star=1, # Lsun
                      R_star=1,
                      heat='radiation') # radiation/accretion
model.get_gasdensitycontrol(abundance=1e-10, # abundance of CH3OH compared to H2
                            snowline=100, # snowline temperature
                            enhancement=1e5, # enhancement factor of abundance inside snowline
                            gas_inside_rcb=True)

##############################################

condition_parms = general_parameters(
    nodust      = False,
    scat        = True,
    extract_gas = True,
)


simulate_mutual_parms = {
    "incl"      : 73,
    "line"      : 240,
    "npix"      : 500,
    "sizeau"    : 200,
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

# class plot_radmc:
#     def __init__(self, **kwargs):
#         for k, v in kwargs.items():
#           # add parameters as attributes of this object
#           setattr(self, k, v)

#     def __del__(self):
#       pass

#     def add_attributes(self, **kwargs):
#       '''
#       Use this function to set the values of the attributs n1, n2, n3,
#       which are number of pixels in the first, second, and third axes. 
#       '''
#       for k, v in kwargs.items():
#         # add parameters as attributes of this object
#         setattr(self, k, v)
    
#     def read_img(self, fname=None,):
#       return
       
