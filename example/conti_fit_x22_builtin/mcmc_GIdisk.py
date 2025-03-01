import numpy as np
from matplotlib import pyplot as plt
import emcee
import corner
import sys
sys.path.insert(0,'../../')
from X22_model.disk_model import *
from radmc.setup import *
from multiprocessing import Pool
import time 
from datetime import datetime
import shutil
import os
from find_center import find_center
import astropy.constants as const
au = const.au.cgs.value
Msun = const.M_sun.cgs.value
yr = 365*24*3600

cb68_alma_list = [
    # '/run/media/hyp0515/fd14f880-ba6f-450f-b82d-98ba3710dc5f/backup/CB68/cb68_alma/cb68_setup1_all.rob2.I.image.tt0.pbcor.smooth.dropdeg.fits',
    # '/run/media/hyp0515/fd14f880-ba6f-450f-b82d-98ba3710dc5f/backup/CB68/cb68_alma/cb68_setup2_all.rob2.I.image.tt0.pbcor.smooth.dropdeg.fits',
    # '/run/media/hyp0515/fd14f880-ba6f-450f-b82d-98ba3710dc5f/backup/CB68/cb68_alma/cb68_setup3_all.rob2.I.image.tt0.pbcor.smooth.dropdeg.fits',
    '/run/media/hyp0515/fd14f880-ba6f-450f-b82d-98ba3710dc5f/backup/CB68_eDisk/CB68_SBLB_continuum_robust_0.0.image.tt0.fits'
]

lambda_list = [
    # 0.13,
    # 0.12,
    # 0.32,
    0.13,
]

sigma_list = [
    # 30e-6,
    # 40e-6,
    # 26e-6,
    21e-6
]


images = []
pa = []
for i, fname in enumerate(cb68_alma_list):
    ra_deg, dec_deg, disk_pa = find_center(fname)
    DI_alma = DiskImage(
        fname = fname,
        ra_deg = ra_deg,
        dec_deg = dec_deg,
        distance_pc = 140,
        rms_Jy = sigma_list[i], # convert to Jy/beam
        disk_pa = disk_pa,
        img_size_au = 60,
        remove_background=True
    )
    images.append(DI_alma)
    pa.append(disk_pa) 



n_processes = 10
nwalkers = 20  # Total number of walkers
ndim = 3        # Dimension of parameter space
niter = 100000     # Number of iterations

def GIdisk(theta):
    # temp_dir_name = f'./temp/{datetime.now().strftime("%H%M%S")}_{datetime.now().microsecond}_{os.getpid()}'
    # os.makedirs('./temp/', exist_ok=True)
    # os.makedirs(temp_dir_name)
    # os.chdir(temp_dir_name)
    
    amax, mstar, mdot = theta
    
    opacity_table = generate_opacity_table_x22(
        a_min=1e-6, a_max=10**(amax-1), # min/max grain size
        q=-3.5, # slope for dust size distribution, dn/da ~ a^q
        dust_to_gas=0.01 # dust-to-gas ratio before sublimation
    )
    # print('pass1')
    disk_property_table = generate_disk_property_table(opacity_table = opacity_table)
    
    D = DiskFitting('CB68', opacity_table, disk_property_table)
    
    cosI = np.cos(np.deg2rad(73))
    D.set_cosI(cosI=cosI)    

    for i, image_class in enumerate(images):
        D.add_observation(image_class, lambda_list[i])    
    ll = -D.mcmc_fit(Mstar=mstar*Msun, Mdot=10**(mdot)*(Msun/yr))
    
    # os.chdir("../..")
    # shutil.rmtree(temp_dir_name)
    
    return ll

# def log_likelihood(theta):
#     return GIdisk(theta=theta)

def log_prior(theta):
    amax, mstar, mdot = theta
    if -3 < amax < 1.5 and 0.1 < mstar < 0.2 and np.log10(2e-7) < mdot < np.log10(8e-7):
        return 0.0
    return -np.inf

def log_probability(theta):
    lp = log_prior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + GIdisk(theta=theta)


pos = [np.array([-0.5, 0.15, np.log10(5e-7)]) + [1e-1, 1e-3, 1e-2] * np.random.randn(ndim) for i in range(nwalkers)]

# File for saving progress
progress_file = "progress.h5"
backend = emcee.backends.HDFBackend(progress_file)
backend.reset(nwalkers, ndim)

# Initialize the sampler
with Pool(n_processes) as pool:
    sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability, pool=pool, backend=backend, a=2)
    sampler.run_mcmc(pos, niter, progress=True)