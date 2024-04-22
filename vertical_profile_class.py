import numpy as np
from disk_model import *
from scipy.interpolate import interp2d
# All unit are in cgs except height and radius are in au
#########################################################################################

class DiskModel_vertical:

    def __init__(self, opacity_table, disk_property_table):
        self.DM_horizontal = DiskModel(opacity_table, disk_property_table)
        return
    
    def input_disk_parameter(self, Mstar, Mdot, Rd, Q, N_R, cut_r_min=True):
        """
        Run Wenrui's radial profile to get initial r-dependent profiles.

        Args:
            Mstar     : mass of protostar
            Mdot      : rate of mass infall from the envelope onto the disk
            Rd        : radius of the disk
            Q         : Toomre index
            N_R       : resolution of radius grid
            cut_r_min : cut_r_min is to cut the data inside the R_min where Sigma = 0 which makes some calculation errors
        """
        self.Mstar = Mstar
        self.Mdot = Mdot
        self.Rd = Rd        
        self.NR = N_R        
        self.DM_horizontal.generate_disk_profile(Mstar=self.Mstar, Mdot=self.Mdot, Rd=self.Rd, Q=Q, N_R=self.NR)
        self.R_grid = self.DM_horizontal.R[1:]/au  # radial grid
        self.T_mid = self.DM_horizontal.T_mid  # temperature at midplane(z=0)
        self.T_eff = self.DM_horizontal.T_eff  # effective temperature extracted from Wenrui's disk_model
        self.M = self.DM_horizontal.Sigma/2  # Sigma/2 is the definition of M in Hubeny+90
        self.tau_r_mid = self.DM_horizontal.tau_r_mid
        self.tau_p_mid = self.DM_horizontal.tau_p_mid
        if cut_r_min:
            cut = np.argmax(self.M > 0)
            self.R_grid = self.R_grid[cut:]
            self.T_mid = self.T_mid[cut:]
            self.T_eff = self.T_eff[cut:]
            self.M = self.M[cut:]
            self.tau_p_mid = self.tau_p_mid[cut:]
            self.tau_r_mid = self.tau_r_mid[cut:]
            self.NR = len(self.R_grid)
            self.cut = cut
        self.Q = G*Mstar*(self.R_grid*au)**(-3)  # (2.2) effective vertical gravity
        self.make_position_map()
        return
    
    def make_position_map(self):
        """
        Make vertical grid
        """
        self.NZ = self.NR
        Z_max = self.Rd
        Z_grid = np.append(np.logspace(np.log10(Z_max/au), np.log10(0.001), self.NZ-1), 0)
        self.Z_grid = Z_grid[::-1]
        R, Z = np.meshgrid(self.R_grid, self.Z_grid, indexing = 'ij')
        pos_map = np.dstack((R, Z))
        self.pos_map = pos_map
        self.precompute_property()
        return
    
    def precompute_property(self, miu=2.3, factor=1):
        """
        To calculate the gas pressure scale height and the radiation pressure scale height
        """
        def cg(miu, factor,  T):  # (5.4) sound speed assiciated with the gas pressure
            """
            miu   : mean molecular weight
            factor: N/(N-n_e)
            T     : temperature

            "In the initialization step, c_g is determined by eq(5.4) with T = T_eff, 
            and with some suitable estimate for N/(N-n_e)"    
            """
            return (kB*T*factor/(miu*mp))**0.5       

        def H_g(cg, Q):  # (4.4) gas pressure scale height
            H_g = (2*cg**2/Q)**0.5
            return H_g/au        

        def H_r(T, Q):  # (4.5) radiation pressure scale height
            get_kappa = self.DM_horizontal.get_kappa_r
            kappa_r = get_kappa(T)
            H_r = (sigma_SB/c_light)*(kappa_r/Q)*T**4
            return H_r/au, kappa_r

        self.cg = cg(miu, factor, self.T_eff)
        self.H_g = H_g(self.cg, self.Q)
        self.H_r, self.kappa_r = H_r(self.T_eff, self.Q)
        self.make_rho_and_m_map()
        return
        
    def make_rho_and_m_map(self):
        """
        Calculate volumn density and mass-depth scale
        """
        M = self.M*1
        h_grid = np.empty((self.NR))
        rho_map = np.empty((self.NR, self.NZ))
        m_map = np.empty((self.NR, self.NZ))
        for r in range(self.NR):
            z_grid = self.Z_grid
            h_g = self.H_g[r]
            h_r = self.H_r[r]
            """
            A simple method to find h: 
            Different h can cause the last element of the m_grid unequal to M, 
            so h can be set as an array from 0 to any large number, iterate until the last element of the m_grid
            nearly equal to M, check which h makes the closest value, and set the value to be the h at that radius.
            """
            h_itr = np.linspace(0.01, 100, 1000000)
            err = 1
            h_index = 0

            while err > 0.1:
                h = h_itr[h_index]
                rho_0 = M[r]/h
                rho_grid = np.where(z_grid <=h,  # condition
                                    rho_0*np.exp(-(z_grid**2/h_g**2)*(1-(h_r/h))),  # if true
                                    rho_0*np.exp(-((z_grid-h_r)/h_g)**2)*np.exp(-((h-h_r)/h_g)*h_r/h_g))  # if false
                rho_grid = np.maximum(rho_grid, 1e-10)  # set lower limit of density
                dz = z_grid*1
                dz[1:] = z_grid[1:]-z_grid[:-1]
                m_grid = np.cumsum(rho_grid*dz)  # (2.3) # m_1 sufficiently small and m_ND = M
                err = np.abs(1-m_grid[-1]/M[r])
                h_index +=1                
            h_grid[r] = h_itr[h_index-1]
            rho_map[r, :] = rho_grid
            m_map[r, :] = m_grid[::-1]

        self.rho_map = rho_map/au  # the divided by au is the made-up of units errors
        self.m_map = m_map
        self.H = h_grid
        self.make_tau_and_T_map()
        return

    
    def make_tau_and_T_map(self):
        """
        Using disk_property_table to calculate kappa_r from T, and further calculating tau_r and T
        """
        
        def get_kappa_from_T(T):  # 2(d)
            # using the functions below can avoid using saha equation to calculate n_e, the step described on 2(d)
            get_kappa_r = self.DM_horizontal.get_kappa_r  # interpolating function extracted from Wenrui's Disk_Model class
            kappa_r = get_kappa_r(T)
            return kappa_r
        def get_tau_from_kappa(m_grid, kappa_r):  # 2(a)
            dm = m_grid.copy()
            dm[1:] = m_grid[1:]-m_grid[:-1]
            dm = dm[:len(kappa_r)]
            dtau_r = dm*kappa_r  # (3.4)
            # The upper and lower bound are m and 0 in m-variable,
            # meaning that the upper and lower bound in z-variable are z and infinity.
            tau_r = np.cumsum(dtau_r)[-1]
            return tau_r
        def get_T_from_tau(r, tau, T_eff):  # 2(b)
            tau_p = self.tau_p_mid[r]*2
            tau_r = self.tau_r_mid[r]*2
            T = (T_eff**4*(3/4)*(tau*(1-tau/tau_r)+(1/np.sqrt(3))+(1/(1.5*tau_p))))**(1/4)  # (3) from Wenrui+22
            return T
        
        T_map = np.empty((self.NR, self.NZ))
        tau_r_map = np.empty((self.NR, self.NZ))
        kappa_r_map = np.empty((self.NR, self.NZ))        
        for r in range(self.NR):
            t_eff = self.T_eff[r]
            m_grid = self.m_map[r, ::-1]
            kappa_r_grid = np.array([])
            T_grid = np.array([])
            tau_r_grid = np.array([])
            for z in range(self.NZ):

                if z == 0:  # initialize kappa_r_grid
                    T_old = 5  # the lowest temperature in Wenrui's code
                    kappa_r_grid = np.append(kappa_r_grid, get_kappa_from_T(T_old))
                else:
                    T_old = T_grid[-1]

                err = 1
                n_itr = 0
                n_itr_max = 50
                # iterate until T_new = T_old
                while err > 1e-10 and n_itr < n_itr_max:
                    kappa_r_grid[z] = get_kappa_from_T(T_old)
                    tau_r = get_tau_from_kappa(m_grid, kappa_r_grid)
                    T_new = get_T_from_tau(r, tau_r, t_eff)
                    err = np.abs(1-T_new/T_old)
                    n_itr +=1
                    T_old = T_new
                T_grid = np.append(T_grid, T_new)
                tau_r_grid = np.append(tau_r_grid, tau_r)
                kappa_r_grid = np.append(kappa_r_grid, get_kappa_from_T(T_new))
            T_map[r, :] = T_grid[::-1]
            tau_r_map[r, :] = tau_r_grid[::-1]
            kappa_r_grid = kappa_r_grid[1:]
            kappa_r_map[r, :] = kappa_r_grid[::-1]
        
        self.T_map = T_map
        self.tau_r_map = tau_r_map
        self.kappa_r_map = kappa_r_map
        return

    def pancake_model(self):
        rho_pancake = 1e-18*np.ones((self.NR, 10))
        T_pancake = 100*np.ones((self.NR, 10))
        self.rho_map = np.append(rho_pancake, np.zeros((self.NR, self.NZ-10)), axis=1)
        self.T_map = np.append(T_pancake, np.zeros((self.NR, self.NZ-10)), axis=1)
        return

    def extend_to_spherical(self, NTheta, NPhi):
        self.NTheta = NTheta
        self.NPhi = NPhi  # Since this model is only axisymmetric so the value of NPhi isn't important

        pos_map = self.pos_map.copy()        
        r_map = np.sqrt(pos_map[:, :, 0]**2+ pos_map[:, :, 1]**2) # distance map of every points
        r_min = np.min(r_map)
        self.rmin = r_min
        r_grid = np.logspace(np.log10(r_min), np.log10(self.Rd/au), self.NR)
        self.r_sph = r_grid

        theta_map = np.arccos(pos_map[:, :, 1]/r_map)
        theta_min = np.deg2rad(30) # the starting angle of theta
        theta_grid = np.logspace(np.log10(theta_min), np.log10(np.max(theta_map)), NTheta)
        theta_grid = -1*theta_grid + 0.5*np.pi + theta_min
        theta_grid = theta_grid[::-1]       
        theta_grid_down = -theta_grid[::-1] + np.pi
        self.theta_sph = np.concatenate((theta_grid, theta_grid_down))     
    
        """
        using interpolation to project data from cylindrical coordinates to 
        spherical coordinates
        """
        z_sph_in_cyl = r_grid[:, np.newaxis] * np.cos(theta_grid) # decreasing when theta increases
        r_sph_in_cyl = r_grid[:, np.newaxis] * np.sin(theta_grid) # increasing when theta increases
        '''
        sph_in_cyl[i, j] : 2D array
        i                : index of r_sph
        j                : index of theta_sph
        '''
        def interpolate(data_map):  
            interpolator = interp2d(self.Z_grid, self.R_grid, data_map, kind='linear')
            interpolated_data = np.empty((self.NR, self.NTheta))
            for r in range(self.NR):
                for theta in range(self.NTheta):
                    interpolated_data[r, theta] = interpolator(z_sph_in_cyl[r, theta], r_sph_in_cyl[r, theta])            
            return interpolated_data    
        def mirror_with_r_plane(map):
            map_mirror = np.fliplr(map)            
            return np.concatenate((map[:, :-1], map_mirror), axis= 1)        
        def rotate_around_theta_axis(map):
            map_3d = np.tile(map[:, :, np.newaxis], (1, 1, self.NPhi))            
            return map_3d        
        rho_sph_2d = interpolate(self.rho_map)
        T_sph_2d = interpolate(self.T_map)        

        """
        Mask the section where the limitation of Wenrui's code encounters
        """
        mask_condition = r_sph_in_cyl < r_min
        rho_sph_2d[mask_condition] = 1e-5/au
        T_sph_2d[mask_condition] = 1200

        self.rho_sph = rotate_around_theta_axis(mirror_with_r_plane(rho_sph_2d)) 
        self.T_sph = rotate_around_theta_axis(mirror_with_r_plane(T_sph_2d))
        self.make_spherical_grid()
        return

    def make_spherical_grid(self):
        """
        Making the grid which RADMC3D adopts
        """
        def make_boundary(grid):
            grid = grid.copy()
            difference = (grid[1:] - grid[:-1])/2
            difference = np.append(difference, difference[-1])
            grid = grid - difference
            return np.append(grid, grid[-1] + 2*difference[-1])
        self.r_sph_grid = make_boundary(self.r_sph)*au
        theta_sph_grid = make_boundary(self.theta_sph)
        theta_sph_grid = np.delete(theta_sph_grid, self.NTheta)
        self.theta_sph_grid = theta_sph_grid
        self.phi_sph_grid = np.linspace(0, 2*np.pi, self.NPhi+1)        
        return