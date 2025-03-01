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
import h5py
from emcee.moves import StretchMove, GaussianMove


n_processes = 20
nwalkers = 20  # Total number of walkers
ndim = 3        # Dimension of parameter space
niter = 100000     # Number of iterations

# # CB68's observation
observed_Freq = np.array([
                        # 95.0,
                        233.8,
                        246.7
                        ])
observed_Flux = np.array([
                        # 7.6,
                        59, 
                        66
                        ])
err           = np.array([
                        # 0.022,
                        0.029,
                        0.040
                        ])



def sed_model(theta):
    temp_dir_name = f'./temp/{datetime.now().strftime("%H%M%S")}_{datetime.now().microsecond}_{os.getpid()}'
    os.makedirs('./temp/', exist_ok=True)
    os.makedirs(temp_dir_name)
    os.chdir(temp_dir_name)
    
    #lam3 = cc*1e4/(95.00*1e9)
    lam2 = cc*1e4/(233.8*1e9)
    lam1 = cc*1e4/(246.7*1e9)
    # lam12    = np.logspace(np.log10(lam1),np.log10(lam2),2,endpoint=False)
    # lam23    = np.logspace(np.log10(lam2),np.log10(lam3),3,endpoint=True)
    # lam      = np.concatenate([lam1,lam2])
    lam      = np.logspace(np.log10(lam1),np.log10(lam2),2,endpoint=True)
    nlam     = lam.size
    with open('camera_wavelength_micron.inp', 'w+') as f:
        f.write('%d\n'%(nlam))
        for value in lam:
            f.write('%13.6e\n'%(value))
    
    amax, Mstar, Mdot = theta
 
    model = radmc3d_setup(silent=False)
    model.get_mastercontrol(filename=None,
                            comment=None,
                            incl_dust=1,
                            incl_lines=1,
                            nphot=500000,
                            nphot_scat=5000000,
                            nphot_spec=500000, 
                            scattering_mode_max=2,
                            istar_sphere=1,
                            num_cpu=1)
    model.get_linecontrol(filename=None,
                        methanol='ch3oh leiden 0 0 0')
    model.get_continuumlambda(filename=None,
                            comment=None,
                            lambda_micron=None,
                            append=False)
    model.get_diskcontrol(  d_to_g_ratio = 0.01,
                            a_max=10**amax, 
                            Mass_of_star=Mstar, 
                            Accretion_rate=10**Mdot,
                            Radius_of_disk=30,
                            NR=200,
                            NTheta=200,
                            NPhi=10)
    model.get_vfieldcontrol(Kep=True,
                            vinfall=None,
                            Rcb=None,
                            outflow=None)
    model.get_heatcontrol(heat='accretion')
    model.get_gasdensitycontrol(abundance=1e-10,
                                snowline=None,
                                enhancement=1e5,
                                gas_inside_rcb=True)

    os.system(f'radmc3d spectrum incl 73 loadlambda noline noscat')
    s = readSpectrum('spectrum.out')
    lam = s[:, 0]
    nu = (1e-2*cc)*1e-9/(1e-6*lam) # GHz
    fnu = s[:, 1]*1e26/(140**2) # mJy

    os.chdir("../..")
    shutil.rmtree(temp_dir_name)
    
    return nu[::-1], fnu[::-1]

nu_record  = []
fnu_record = [] 

# Log-prob, log-likelihood, and prior
def log_likelihood(theta, Freq, Flux, err):

    freq_model, flux_model = sed_model(theta=theta)
    nu_record.append(freq_model)
    fnu_record.append(flux_model)
    np.savez('record.npz',
         nu  = np.array(nu_record),
         fnu = np.array(fnu_record))
    
    mapping = {95.0: 0, 233.8: 0, 246.7: 1}
    freq_index = [mapping[value] for value in Freq]
    
    return -0.5 * np.sum((Flux - flux_model[freq_index]) ** 2 / err**2)

def log_prior(theta):
    a_max, Mstar, Mdot = theta
    if -3 < a_max < 3 and 0.1 < Mstar < 0.2 and -8 < Mdot < -5.5:
        return 0.0
    return -np.inf

def log_probability(theta, Freq, Flux, err):
    lp = log_prior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood(theta, Freq, Flux, err)


pos = [np.array([-2, 0.14, -6.5]) + [1e-1, 1e-2, 1e-1] * np.random.randn(ndim) for i in range(nwalkers)]

# moves = [
#     (StretchMove(a=2.0), 0.8),  # StretchMove with 80% probability
#     (GaussianMove(cov=0.5 * np.identity(ndim)), 0.2),  # GaussianMove with 20% probability
# ]

# File for saving progress
progress_file = "progress.h5"
backend = emcee.backends.HDFBackend(progress_file)
# backend.reset(nwalkers, ndim)

# Initialize the sampler
with Pool(n_processes) as pool:
    sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability, args=(observed_Freq, observed_Flux, err), pool=pool, backend=backend, a=2)
    # sampler.run_mcmc(pos, niter, progress=True)
    sampler.run_mcmc(None, niter, progress=True)
    

np.savez('record.npz',
         nu  = np.array(nu_record),
         fnu = np.array(fnu_record))

