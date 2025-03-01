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

reader = emcee.backends.HDFBackend("progress.h5")


fig, axes = plt.subplots(3, figsize=(10, 10), sharex=True)
samples = reader.get_chain()
label = [r'log($a_{max}$)',r'$M_{*}$', r'log($\dot{M}$)']
for i in range(len(label)):
    ax = axes[i]
    ax.plot(samples[:, :, i], "k", alpha=0.3)
    ax.set_xlim(0, len(samples))
    ax.set_ylabel(label[i])
    ax.yaxis.set_label_coords(-0.1, 0.5)
# plt.show()
plt.savefig('chain_step.pdf', transparent = True)
plt.close()



# flat_samples = reader.get_chain(flat=True, discard=100)
chains = reader.get_chain(discard=500, flat=False)
n_steps, n_walkers, n_params = chains.shape

# Calculate acceptance fractions manually
# Count accepted steps for each walker
acceptance_fractions = np.mean(
    np.diff(reader.get_log_prob(discard=500, flat=False), axis=0) != 0, axis=0
)

# Define a threshold for stuck walkers
threshold = 0.1
valid_walkers = np.where(acceptance_fractions > threshold)[0]

# Filter the chains to exclude stuck walkers
filtered_chains = chains[:, valid_walkers, :]

# Flatten the filtered chains
flat_samples = filtered_chains.reshape(-1, n_params)[::1]
fig = corner.corner(
    flat_samples[::], labels=label,
    show_titles=True, plot_datapoints=True, quantiles=[0.16, 0.5, 0.84])
fig.savefig('posterior.pdf', transparent=True)
plt.close()



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

model = np.load('record.npz')
nu = model['nu'][::]
fnu = model['fnu'][::]


for i in range(len(nu)):

    plt.plot(nu[i], fnu[i], "C1", alpha=0.01)
plt.errorbar(observed_Freq, observed_Flux, yerr=err, fmt=".k", capsize=0)
plt.xlabel('$\\nu [GHz]$')
plt.ylabel('Flux density [mJy]')
plt.ylim((50, 70))
plt.savefig('plot.pdf', transparent=True)
plt.close()
