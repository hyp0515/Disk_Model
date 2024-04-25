import numpy as np
from mpl_toolkits.mplot3d import axes3d
from matplotlib import pyplot as plt
#
# Some natural constants
#
au  = 1.49598e13     # Astronomical Unit       [cm]
pc  = 3.08572e18     # Parsec                  [cm]
ms  = 1.98892e33     # Solar mass              [g]
ts  = 5.78e3         # Solar temperature       [K]
ls  = 3.8525e33      # Solar luminosity        [erg/s]
rs  = 6.96e10        # Solar radius            [cm]
#
# Monte Carlo parameters
#
nphot    = 1000000
#
# Disk Model
#
from disk_model import *
from vertical_profile_class import DiskModel_vertical
class problem_setup:
    def __init__(self, a_max, Mass_of_star, Accretion_rate, Radius_of_disk, pancake=False):        
        opacity_table = generate_opacity_table(a_min=0, a_max=a_max, q=-3.5, dust_to_gas=0.01)
        disk_property_table = generate_disk_property_table(opacity_table)

        DM = DiskModel_vertical(opacity_table, disk_property_table)
        DM.input_disk_parameter(Mstar=Mass_of_star, Mdot=Accretion_rate,
                                Rd=Radius_of_disk, Q=1.5, N_R=100)
        if pancake is True:
            DM.pancake_model()
        DM.extend_to_spherical(NTheta=200, NPhi=200)
        self.r_sph = DM.r_sph
        self.theta_sph = np.delete(DM.theta_sph, DM.NTheta)
        #
        # Write the wavelength_micron.inp file
        #
        # lam1     = 0.1e0
        # lam2     = 7.0e0
        # lam3     = 25.e0
        # lam4     = 1.0e4
        lam1     = 3.0e2
        lam2     = 3.0e3
        lam3     = 3.0e4
        lam4     = 3.0e5
        n12      = 20
        n23      = 100
        n34      = 30
        lam12    = np.logspace(np.log10(lam1),np.log10(lam2),n12,endpoint=False)
        lam23    = np.logspace(np.log10(lam2),np.log10(lam3),n23,endpoint=False)
        lam34    = np.logspace(np.log10(lam3),np.log10(lam4),n34,endpoint=True)
        lam      = np.concatenate([lam12,lam23,lam34])
        nlam     = lam.size
        with open('wavelength_micron.inp', 'w+') as f:
            f.write('%d\n'%(nlam))
            for value in lam:
                f.write('%13.6e\n'%(value))
        #
        # Star parameters and Write the star.inp file
        #
        mstar    = DM.Mstar
        rstar    = rs
        tstar    = ts
        pstar    = np.array([0.,0.,0.])
        with open('stars.inp','w+') as f:
            f.write('2\n')
            f.write('1 %d\n\n'%(nlam))
            f.write('%13.6e %13.6e %13.6e %13.6e %13.6e\n\n'%(rstar,mstar,pstar[0],pstar[1],pstar[2]))
            for value in lam:
                f.write('%13.6e\n'%(value))
            f.write('\n%13.6e\n'%(-tstar))
        #
        # Write the grid file
        #
        iformat     = 1
        grid_style  = 0
        coordsystem = 150  # 100 <= coordsystem < 200 is spherical
        grid_info   = 0  # advised to set =0
        incl_r      = 1
        incl_theta  = 1
        incl_phi    = 1
        nr          = DM.NR
        ntheta      = 2*DM.NTheta-1 # NTheta is semispherical
        nphi        = DM.NPhi
        with open('amr_grid.inp', "w") as f:
            f.write(str(iformat)+'\n')
            f.write(str(grid_style)+'\n')
            f.write(str(coordsystem)+'\n')
            f.write(str(grid_info)+'\n')
            f.write('%d %d %d\n'%(incl_r, incl_theta, incl_phi))
            f.write('%d %d %d\n'%(nr, ntheta, nphi))
            for value in DM.r_sph_grid:
                f.write('%13.13e\n'%(value))
            for value in DM.theta_sph_grid:
                f.write('%13.13e\n'%(value))
            for value in DM.phi_sph_grid:
                f.write('%13.13e\n'%(value))
        #
        # Write the density file
        #
        nspec       = 1
        with open('dust_density.inp', "w") as f:
            f.write(str(iformat)+'\n')
            f.write('%d\n'%(nr*ntheta*nphi))
            f.write(str(nspec)+'\n')
            data = 0.01*DM.rho_sph.ravel(order='F')         # Create a 1-D view, fortran-style indexing
            data.tofile(f, sep='\n', format="%13.6e")
            f.write('\n')
        #
        # Write the temperature file
        #
        with open('dust_temperature.dat', "w") as f:
            f.write(str(iformat)+'\n')
            f.write('%d\n'%(nr*ntheta*nphi))
            f.write(str(nspec)+'\n')
            data = DM.T_sph.ravel(order='F')         # Create a 1-D view, fortran-style indexing
            data.tofile(f, sep='\n', format="%13.6e")
            f.write('\n')
        #
        # Write the radmc3d.inp control file
        #
        with open('radmc3d.inp','w+') as f:
            f.write('nphot = %d\n'%(nphot))
            f.write('scattering_mode_max = 0\n')   # Put this to 1 for isotropic scattering
            # f.write('iranfreqmode = 1\n')
            f.write('istar_sphere = 1\n')
            f.write('tgas_eq_tdust = 1\n')
            f.write('setthreads = 8\n')
        #
        # Write dust opacity files
        #
        iformat   = 3
        nlam      = len(opacity_table['lam'])
        lam       = opacity_table['lam']*1e4     # lam in opacity_table is in cgs while RADMC3D uses micro
        kappa_abs = opacity_table['kappa']
        kappa_sca = opacity_table['kappa_s']
        g         = opacity_table['g']
        # for idx, composition in enumerate(['water','silicate','troilite','refractory_organics']):
        for idx, composition in enumerate(['silicate']): # for now, only silicate is considered
            with open('dustkappa_'+composition+'.inp', "w+") as f:
                f.write(str(iformat)+'\n')
                f.write(str(nlam)+'\n')
                for lam_idx in range(nlam):
                    f.write('%13.6e %13.6e %13.6e %13.6e\n'%(lam[lam_idx],kappa_abs[idx,lam_idx],kappa_sca[idx,lam_idx],g[idx,lam_idx]))
        #
        # Dust opacity control file
        #
        with open('dustopac.inp','w+') as f:
            f.write('2                          Format number of this file\n')
            f.write('1                          Nr of dust species\n')
            f.write('============================================================================\n')
            f.write('1                          Way in which this dust species is read\n')
            f.write('0                          0=Thermal grain\n')
            f.write('silicate                   Extension of name of dustkappa_***.inp file\n')
        return
