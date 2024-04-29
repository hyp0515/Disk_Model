import numpy as np
from matplotlib import pyplot as plt
from radmc3dPy.image import *
from radmc3dPy.analyze import *
from disk_model import *
from vertical_profile_class import DiskModel_vertical
from problem_setup import problem_setup
from scipy.optimize import curve_fit


###############################################################################
'''
Checking the consistency between X22 model and its extension
'''
# Density (cylindrical coordinates)
# p = problem_setup(a_max=0.01, Mass_of_star=0.14*Msun, Accretion_rate=0.14e-5*Msun/yr, Radius_of_disk=30*au, pancake=False)
# r = p.DM.R_grid
# plt.plot(r, p.DM.M)
# plt.plot(r, p.DM.m_map[:, 0])
# plt.legend(['M', 'm'])
# plt.yscale('log')
# plt.show()
# plt.close()
# z_grid = p.DM.Z_grid
# dz = np.diff(z_grid, prepend=z_grid[0])
# c = np.sum(p.DM.rho_map[20, :]*dz*au)
# print(c, p.DM.m_map[20, 0], p.DM.M[20])



# Temperature (cylindrical coordinates)
# p = problem_setup(a_max=0.01, Mass_of_star=0.14*Msun, Accretion_rate=0.14e-5*Msun/yr, Radius_of_disk=30*au, pancake=False)
# r = p.DM.R_grid
# print(p.DM.T_mid[10], p.DM.T_map[10, 0])
# plt.plot(r, p.DM.T_mid)
# plt.plot(r, p.DM.T_map[:, 0])
# plt.legend(['Tmid', 'Tmap'])
# plt.show()
# plt.close()



# Density (spherical coordinates)
p = problem_setup(a_max=0.01, Mass_of_star=0.14*Msun, Accretion_rate=0.14e-5*Msun/yr, Radius_of_disk=30*au, pancake=False)
rho = p.DM.rho_sph[:, :, 0]

r = p.DM.r_sph*au
dr = np.diff(r, prepend=r[0])
theta = p.DM.theta_sph
dtheta = np.diff(theta, prepend=theta[0])

r_2d = np.tile(r[:, np.newaxis], (1, len(theta)))
dr_2d = np.tile(dr[:, np.newaxis], (1, len(theta)))
theta_2d = np.tile(theta[np.newaxis, :], (len(r), 1))
dtheta_2d = np.tile(dtheta[np.newaxis, :], (len(r), 1))

integrand_rho = rho*r_2d*dr_2d*dtheta_2d
total_rho = np.sum(integrand_rho)

sigma = 2*p.DM.M
R = p.DM.R_grid*au
dR = np.diff(R, prepend=R[0])
integrand_sigma = sigma*dR
total_sigma = np.sum(integrand_sigma)

# print(total_rho/Msun, total_sigma/Msun)

phi = np.linspace(0, 2*pi, 100, endpoint=True)
dphi = np.diff(phi, prepend=phi[0])
dphi_3d = np.tile(dphi[np.newaxis, np.newaxis, :], (len(r), len(theta), 1))

rho_3d = np.tile(rho[:, :, np.newaxis], (1, 1, len(phi)))
r_3d = np.tile(r_2d[:, :, np.newaxis], (1, 1, len(phi)))
dr_3d = np.tile(dr_2d[:, :, np.newaxis], (1, 1, len(phi)))
dtheta_3d = np.tile(dtheta_2d[:, :, np.newaxis], (1, 1, len(phi)))
sintheta = np.sin(theta_2d)
sintheta_3d = np.tile(sintheta[:, :, np.newaxis], (1, 1, len(phi)))

integrand_rho = rho_3d*r_3d*r_3d*sintheta_3d*dr_3d*dtheta_3d*dphi_3d
total_rho_3d = np.sum(integrand_rho)

sigma_2d = np.tile(sigma[:, np.newaxis], (1, len(phi)))
R_2d = np.tile(R[:, np.newaxis], (1, len(phi)))
dR_2d = np.tile(dR[:, np.newaxis], (1, len(phi)))
dphi_2d = np.tile(dphi[np.newaxis, :], (len(R), 1))

integrand_sigma = sigma_2d*dR_2d*R_2d*dphi_2d
total_sigma_3d = np.sum(integrand_sigma)

print(total_rho_3d/Msun, total_sigma_3d/Msun)