import numpy as np
from matplotlib import pyplot as plt
from disk_model import *
from vertical_profile_class import DiskModel_vertical

opacity_table = generate_opacity_table(a_min=0, a_max=0.1, q=-3.5, dust_to_gas=0.01)
disk_property_table = generate_disk_property_table(opacity_table)
DM = DiskModel(opacity_table=opacity_table, disk_property_table=disk_property_table)
DM.generate_disk_profile(Mstar=15*Msun, Mdot=15*1e-5*Msun/yr,Rd=50*au, Q=1.5, N_R=200)
plt.plot(DM.R[1:]/au, DM.Sigma)
DM.generate_disk_profile(Mstar=10*Msun, Mdot=10*1e-5*Msun/yr,Rd=50*au, Q=1.5, N_R=200)
plt.plot(DM.R[1:]/au, DM.Sigma)
DM.generate_disk_profile(Mstar=5*Msun, Mdot=5*1e-5*Msun/yr,Rd=50*au, Q=1.5, N_R=200)
plt.plot(DM.R[1:]/au, DM.Sigma)
plt.show()
# DM = DiskModel_vertical(opacity_table, disk_property_table, Mstar=0.5*Msun, Mdot=5e-6*Msun/yr, Rd=50*au, Z_max=50*au, Q=1.5, N_R=250, N_Z=250)
# DM.precompute_property(miu=2, factor=1.5)
# DM.extend_to_spherical(NTheta=500)

# print(DM.theta_grid[-1])
# print(0.5*np.pi)

# print(DM.r_sph) 