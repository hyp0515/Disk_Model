import numpy as np
import matplotlib.pyplot as plt
from generate_cube import *

from CB68.data_dict import data_dict

methanol_data = data_dict["ch3oh_218_faust"]
beam_axis = [methanol_data["bmaj"], methanol_data["bmin"]]
beam_pa = methanol_data["bpa"]

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
                generate_model( amax     =  a, # mm
                                mstar    =  0.14, # Msun
                                mdot     =  mdot, # Msun/yr
                                Q        =  Q, # Toomre Q
                                snowline =  snowline, # temperature of sublimation
                                rcb      =  None)

                produce_cube(fname=f"a_{a}_mdot_{mdot}_Q_{Q}_snowline_{snowline}",channel=False, pv=True)

# im_conti, im_conv_conti = initialize_image(
#     fname='./test/outfile/pv_test_conti.out',
#     convolve=True)

# im_line, im_conv_line = initialize_image(
#     fname='./test/outfile/pv_test_scat.out',
#     convolve=True)


