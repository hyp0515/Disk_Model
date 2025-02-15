import io
import os
import contextlib
from datetime import datetime
import shutil
import numpy as np

from matplotlib import pyplot as plt
from matplotlib.patches import Ellipse
import emcee
from multiprocessing import Pool

from radmc3dPy import image

import sys
sys.path.insert(0,'../../')
from X22_model.disk_model import *
from radmc.setup import *
from CB68.data_dict import data_dict

from astropy.coordinates import SkyCoord

n_processes = 20
nwalkers = 6  # Total number of walkers
ndim = 3        # Dimension of parameter space
niter = 100000     # Number of iterations

"""
Initialize continuum data
"""
# dust continuum data
conti_data    = [data_dict["1.3_edisk"]]
conti_lam_list = [conti_data[i]["wav"] for i in range(len(conti_data))]
conti_sigma_list  = [conti_data[i]["sigma"] for i in range(len(conti_data))]
disk_posang = 45
desire_size = [50]

centroid_position = [
    ['16h57m19.6429s', '-16d09m24.027s'],
]

observation_data = []
beam_pa = []
beam_axis = []
size_au = []
npix = []

for i, data in enumerate(conti_data):

    # Convert the center position to ra, dec in degree
    c = SkyCoord(centroid_position[i][0], centroid_position[i][1], frame='icrs')
    ra_deg = c.ra.deg
    dec_deg = c.dec.deg

    # Initialize the DiskImage class
    image_class = DiskImage(
        fname = data["fname"],
        ra_deg = ra_deg,
        dec_deg = dec_deg,
        distance_pc = 140,
        rms_Jy = data["sigma"],
        disk_pa = disk_posang,
        img_size_au = desire_size[i],
        remove_background=False,
    )

    observation_data.append(image_class.img)
    beam_pa.append(image_class.beam_pa)
    beam_axis.append([  image_class.beam_maj_au/image_class.distance_pc,
                        image_class.beam_min_au/image_class.distance_pc])
    size_au.append(image_class.img_size_au*2)
    npix.append(image_class.img.shape[0])



# sed data

sed_data = [data_dict["3.2_faust"], data_dict["1.3_faust"], data_dict["1.2_faust"]]
sed_freq_list = [sed_data[i]["freq"] for i in range(len(sed_data))]
sed_flux_list = [sed_data[i]["flux"] for i in range(len(sed_data))]
sed_err_list = [sed_data[i]["sigma"] for i in range(len(sed_data))]


"""
Disk model + radmc3d
---------------------------------------------
This function produce the model image using radmc3d.
RADMC-3D is works in each directory named by the current time and pid
in order to do parallelization.
---------------------------------------------
"""
def disk_model(parms):
    amax, Mdot, Q = parms
    model = radmc3d_setup(silent=True)
    model.get_mastercontrol(filename=None,
                            comment=None,
                            incl_dust=1,
                            incl_lines=1,
                            nphot=500000,
                            nphot_scat=1000000,
                            nphot_spec=100000, 
                            scattering_mode_max=2,
                            istar_sphere=1,
                            num_cpu=6)
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
    model.get_heatcontrol(heat='accretion')

def radmc_conti():
    
    model_image_list = []
    beam_per_pix = []
    for i, wav in enumerate(conti_lam_list):
        os.system(f'radmc3d image npix {npix[i]} sizeau {size_au[i]} posang {-disk_posang} incl 73 lambda {wav*1000} noline > /dev/null 2>&1')
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            im = image.readImage()
        # Convolve the image with the beam
        model_image = im.imConv(dpc=140, fwhm=beam_axis[i], pa=-beam_pa[i])
        pixel_area = (size_au[i]/npix[i]/140)**2
        beam_area = beam_axis[i][0]*beam_axis[i][1]*np.pi/(4*np.log(2))
        beam_per_pix.append(pixel_area/beam_area)
        model_image_list.append(model_image.imageJyppix[:, :, 0]/(140**2)*(beam_area/pixel_area))
    return model_image_list, beam_per_pix

def radmc_sed():
    lam = np.array([cc*1e4/(freq*1e9) for freq in sed_freq_list])
    nlam     = lam.size
    with open('camera_wavelength_micron.inp', 'w+') as f:
        f.write('%d\n'%(nlam))
        for value in lam:
            f.write('%13.6e\n'%(value))
    os.system(f'radmc3d spectrum incl 73 loadlambda noline')
    s = readSpectrum('spectrum.out')
    # lam = s[:, 0]
    # nu = (1e-2*cc)*1e-9/(1e-6*lam) # GHz
    fnu = s[:, 1]*1e26/(140**2) # mJy
    return sed_freq_list, fnu[::-1]

def conti_model(theta):
    # Create a temporary directory for model computation
    temp_dir_name = f'./temp/{datetime.now().strftime("%H%M%S")}_{datetime.now().microsecond}_{os.getpid()}'
    os.makedirs('./temp/', exist_ok=True)
    os.makedirs(temp_dir_name)
    os.chdir(temp_dir_name)
    
    disk_model(theta)

    model_image_list, beam_per_pix = radmc_conti()
    freq, flux_model = radmc_sed()


    # Clean up temporary directory
    os.chdir("../..")
    shutil.rmtree(temp_dir_name)
    
    return model_image_list, beam_per_pix, freq, flux_model

"""
log_likelihood, log_prior, log_probability
"""
def log_likelihood(theta, conti_observation, sed_observation, conti_err, sed_err):
    # Compute the model image
    model_image, beam_per_pix, freq, flux_model = conti_model(theta=theta)


    log_likelihood = []
    
    """
    This log-likelihood function is modified from disk_model.py
    where chi-square is chosen to be the minimum of the two chi-square
    considering the observation error and the model error.
    """
    def ll_conti(observation, model, beam_per_pix, sigma_obs, sigma_log_model = np.log(2)/2):

        # These chi-square functions are still needed to be confirmed

        chisq1 = (observation-model)**2/(2*sigma_obs**2)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            chisq2 = np.log(observation/model)**2 / (2*sigma_log_model**2)
        chisq2 = np.nan_to_num(chisq2, nan=1e6)
        chisq = np.minimum(chisq1, chisq2)
        dlog_likelihood =  - chisq
        log_likelihood = np.sum(dlog_likelihood*beam_per_pix)
        return log_likelihood

    nu_record  = []
    fnu_record = [] 

    def ll_sed(freq, flux_observe, flux_model, err):
        nu_record.append(freq)
        fnu_record.append(flux_model)
        np.savez('record.npz',
            nu  = np.array(nu_record),
            fnu = np.array(fnu_record))
        return -0.5 * np.sum((np.array(flux_observe - flux_model) ** 2) / err**2)


    # Calculate the log likelihood for each observation
    for i in range(len(conti_observation)):
        mask = conti_observation[i] < 10*conti_err[i]
        conti_observation[i][mask] = 0
        model_image[i][mask.T] = 0
        log_likelihood.append(ll_conti(conti_observation[i], model_image[i].T, beam_per_pix[i], conti_err[i]))

    log_likelihood.append(ll_sed(freq, sed_observation, flux_model, sed_err))

    return np.average(log_likelihood, weights=[0.5, 0.5])

def log_prior(theta):
    amax, Mdot, Q = theta
    """
    These priors are chosen relatively wide.
    amax: 1um < amax < 1m
    Mdot: 1e-8 < Mdot < 1e-6 Msun/yr
    Toomre index: 0.5 < Q < 2.5 (gravitationally unstable to stable)
    """
    # Define the prior ranges for the parameters
    if -3 < amax < 3 and np.log10(1e-8) < Mdot < np.log10(1e-6) and 0.5 < Q < 2.5:
        return 0.0
    return -np.inf

def log_probability(theta, conti_observation, sed_observation, conti_err, sed_err):
    # Calculate the log probability
    lp = log_prior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood(theta, conti_observation, sed_observation, conti_err, sed_err)


"""
MCMC
"""
def mcmc():
    # Initialize the starting positions for the walkers
    pos = [np.array([-1, np.log10(5e-7), 1.5]) + [1e-1, 1e-2, 1e-1] * np.random.randn(ndim) for i in range(nwalkers)]
    # File for saving progress
    progress_file = "progress.h5"
    backend = emcee.backends.HDFBackend(progress_file)
    backend.reset(nwalkers, ndim)
    # Initialize the sampler
    with Pool(n_processes) as pool:
        sampler = emcee.EnsembleSampler(nwalkers,
                                        ndim, 
                                        log_probability, 
                                        args=(observation_data,
                                              conti_sigma_list,
                                              sed_flux_list,
                                              sed_err_list), 
                                        pool=pool, 
                                        backend=backend)
        sampler.run_mcmc(pos, niter, progress=True)

mcmc()