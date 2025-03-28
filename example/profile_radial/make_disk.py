import os
import sys
sys.path.insert(0,'../../')
from radmc.setup import radmc3d_setup

def generate_disk(  amax   = 0.1,
                    mstar  = 0.14,
                    mdot   = 1e-7,
                    rd     = 30,
                    Q      = 1.5,
                    l_star = 0.89,
                    r_star = 1,
                    heat   = "accretion",
                    dir="./disk/"):
    os.makedirs(dir, exist_ok=True)
    os.chdir(dir)
    model = radmc3d_setup(silent=False)
    model.get_mastercontrol(filename=None,
                            comment=None,
                            incl_dust=1,
                            incl_lines=1,
                            nphot=5000000,
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
                            append=False,
                            silent=True)

    model.get_diskcontrol(  d_to_g_ratio = 0.01,
                            a_max=amax, 
                            Mass_of_star=mstar, 
                            Accretion_rate=mdot,
                            Radius_of_disk=rd,
                            NR=200,
                            NTheta=200,
                            NPhi=2,
                            Q=Q)
    model.get_vfieldcontrol(Kep=True,
                            vinfall=0.5,
                            Rcb=None,
                            outflow=None)
    model.get_heatcontrol(L_star=l_star,
                          R_star=r_star,
                          heat=heat)
    model.get_gasdensitycontrol(abundance=1e-10,
                                snowline=100,
                                enhancement=1e5,
                                gas_inside_rcb=True)
    os.chdir("..")
    return model
##############################################
