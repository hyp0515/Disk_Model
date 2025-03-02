import io
import os
import contextlib
from datetime import datetime
import shutil
import numpy as np
from matplotlib import pyplot as plt
import emcee
from multiprocessing import Pool

from radmc3dPy import image

import sys
sys.path.insert(0,'../../')
from X22_model.disk_model import *
from radmc.setup import *
from CB68.data_dict import data_dict

n_processes = 20
nwalkers = 8  # Total number of walkers
ndim = 4        # Dimension of parameter space
niter = 100000     # Number of iterations

"""
Initialize observation data
"""
edisk_radial = np.load("edisk_radial.npz")
sigma_obs = data_dict["1.3_edisk"]["sigma"]
i_r_obs = edisk_radial["i_r"]
beam_axis = edisk_radial["beam_axis"]
beam_pa = edisk_radial["beam_pa"]
npix = len(i_r_obs)

size_au = 80
disk_posang = 45
wav = 1.3
distance_pc = 140

pixel_area = (size_au/npix/distance_pc)**2
beam_area = beam_axis[0]*beam_axis[1]*np.pi/(4*np.log(2))

"""
Disk model + radmc3d
---------------------------------------------
This function produce the model image using radmc3d.
RADMC-3D is works in each directory named by the current time and pid
in order to do parallelization.
---------------------------------------------
"""
def rotate_image(image, posang):
    if isinstance(image, np.ndarray)!=True:
        image.imageJyppix= ndimage.rotate(image.imageJyppix, posang, reshape=False, axes=(1, 0))
        image.imageJyppix = np.nan_to_num(image.imageJyppix, nan=0)
        return image.imageJyppix[:, :, 0]
    else:
        image = ndimage.rotate(image, posang, reshape=False, axes=(1, 0))
        image = np.nan_to_num(image, nan=0)
        return image

def radial_intensity(image_array, center, width):
    radial_profile = np.mean(image_array[:, center-width//2:center+width//2], axis=1)
    return radial_profile

def radmc_conti(theta):
    r_star, amax, Mdot, Q = theta

    
    model = radmc3d_setup(silent=True)
    model.get_mastercontrol(filename=None,
                            comment=None,
                            incl_dust=1,
                            incl_lines=1,
                            nphot=100000,
                            nphot_scat=500000,
                            scattering_mode_max=2,
                            istar_sphere=1,
                            num_cpu=10,
                            modified_random_walk=1)
    model.get_linecontrol(filename=None,
                        methanol='ch3oh leiden 0 0 0')
    model.get_continuumlambda(filename=None,
                            comment=None,
                            lambda_micron=None,
                            append=False,
                            silent=True)
    model.get_diskcontrol(  d_to_g_ratio = 0.01,
                            a_max=10**amax, 
                            Mass_of_star=0.14, 
                            Accretion_rate=10**Mdot,
                            Radius_of_disk=30,
                            Q=Q,
                            NR=200,
                            NTheta=200,
                            NPhi=10)
    model.get_heatcontrol(L_star=0.89,
                        R_star=r_star,
                        heat="irradiation")
    os.system(f'radmc3d image npix {npix} sizeau {size_au} posang {-disk_posang} incl 73 lambda {wav*1000} noline > /dev/null 2>&1')
    with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
        im = image.readImage()
    im_conv = im.imConv(dpc=distance_pc, fwhm=beam_axis, pa=-beam_pa)
    im_conv = rotate_image(im_conv, disk_posang)
    i_r = radial_intensity(im_conv, npix//2, 10)
    return i_r*(beam_area/pixel_area)/(distance_pc**2)

def conti_model(params):
    # Create a temporary directory for model computation
    temp_dir_name = f'./temp/{datetime.now().strftime("%H%M%S")}_{datetime.now().microsecond}_{os.getpid()}'
    os.makedirs('./temp/', exist_ok=True)
    os.makedirs(temp_dir_name)
    os.chdir(temp_dir_name)
    
    i_r_model = radmc_conti(params)
    
    # Clean up temporary directory
    os.chdir("../..")
    shutil.rmtree(temp_dir_name)
    
    return i_r_model

"""
log_likelihood, log_prior, log_probability
"""
def log_likelihood(theta, i_r_obs, sigma_obs):
    # Compute the model image
    i_r_model = conti_model(params=theta)
    chisq1 = (i_r_obs-i_r_model)**2/(2*sigma_obs**2)
    sigma_log_model = np.log(2)/2
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        chisq2 = np.log(i_r_obs/i_r_model)**2 / (2*sigma_log_model**2)
    chisq2 = np.nan_to_num(chisq2, nan=1e6)
    chisq = np.minimum(chisq1, chisq2)
    dlog_likelihood =  - chisq
    log_likelihood = np.sum(dlog_likelihood)

    return log_likelihood


def log_prior(theta):
    r_star, amax, Mdot, Q = theta
    """
    These priors are chosen relatively wide.
    r_star: 0.1 < r_star < 10 Rsun
    amax: 1um < amax < 1m
    Mdot: 1e-8 < Mdot < 1e-6 Msun/yr
    Toomre index: 0.5 < Q < 2.5 (gravitationally unstable to stable)
    """
    # Define the prior ranges for the parameters
    if 0.1 < r_star < 10 and -3 < amax < 3 and -8 < Mdot < -6 and 0.5 < Q < 2.5:
        return 0.0
    return -np.inf

def log_probability(theta, observation, err):
    # Calculate the log probability
    lp = log_prior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood(theta, observation, err)

"""
This is a debugging function to check the whole process
"""
def debugger(theta=(1, -1, np.log10(5e-7), 1.5)):

    i_r_model = conti_model(theta=theta)
    def log_likelihood(i_r_obs, sigma_obs):
        # Compute the model image
        chisq1 = (i_r_obs-i_r_model)**2/(2*sigma_obs**2)
        sigma_log_model = np.log(2)/2
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            chisq2 = np.log(i_r_obs/i_r_model)**2 / (2*sigma_log_model**2)
        chisq2 = np.nan_to_num(chisq2, nan=1e6)
        chisq = np.minimum(chisq1, chisq2)
        dlog_likelihood =  - chisq
        log_likelihood = np.sum(dlog_likelihood*pixel_area/beam_area)

        return log_likelihood
    ll = log_likelihood(i_r_obs, sigma_obs)
    print(f"Log likelihood: {ll}")
    plt.figure(figsize=(8, 4))
    plt.plot(edisk_radial["r_axis"], i_r_obs*1e3, label="Observation", linestyle="--")
    plt.plot(edisk_radial["r_axis"], i_r_model*1e3, label="Model")
    plt.xlim((-(size_au//2), size_au//2))
    plt.ylim(bottom=0)
    plt.xlabel("Offset (AU)")
    plt.ylabel("Intensity (mJy/beam)")
    plt.title("Continuum radial profile along the major axis")
    plt.legend()
    plt.show()
    plt.close("all")
        

"""
MCMC
"""
def mcmc():
    # Initialize the starting positions for the walkers
    pos = [np.array([1, -1, np.log10(5e-7), 1.5]) + [1e-1, 1e-1, 1e-2, 1e-1] * np.random.randn(ndim) for i in range(nwalkers)]
    # File for saving progress
    progress_file = "progress.h5"
    backend = emcee.backends.HDFBackend(progress_file)
    backend.reset(nwalkers, ndim)
    # Initialize the sampler
    with Pool(n_processes) as pool:
        sampler = emcee.EnsembleSampler(nwalkers,
                                        ndim, 
                                        log_probability, 
                                        args=(i_r_obs, sigma_obs), 
                                        pool=pool, 
                                        backend=backend)
        sampler.run_mcmc(pos, niter, progress=True)

# debugger()
mcmc()