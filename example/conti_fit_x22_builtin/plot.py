import numpy as np
from matplotlib import pyplot as plt
import matplotlib as mpl
from matplotlib.colors import ListedColormap
top = mpl.colormaps['Reds_r'].resampled(128)
bottom = mpl.colormaps['Blues'].resampled(128)

newcolors = np.vstack((top(np.linspace(0, 1, 128)),
                       bottom(np.linspace(0, 1, 128))))
residual_cmp = ListedColormap(newcolors, name='RedsBlue')
from matplotlib.patches import Ellipse
import emcee
import corner
import sys
sys.path.insert(0,'../../')
from X22_model.disk_model import *
from find_center import find_center

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


flat_samples = reader.get_chain(flat=True, discard=100)
chains = reader.get_chain(discard=100, flat=False)
n_steps, n_walkers, n_params = chains.shape

# Calculate acceptance fractions manually
# Count accepted steps for each walker
acceptance_fractions = np.mean(
    np.diff(reader.get_log_prob(discard=100, flat=False), axis=0) != 0, axis=0
)

# Define a threshold for stuck walkers
threshold = 0.001
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


cb68_alma_list = [
    '/run/media/hyp0515/fd14f880-ba6f-450f-b82d-98ba3710dc5f/backup/CB68_eDisk/CB68_SBLB_continuum_robust_0.0.image.tt0.fits'
]

lambda_list = [
    0.13,
]

sigma_list = [
    21e-6
]

samples = reader.get_chain(flat=True, discard=100)

theta_max = samples[np.argmax(reader.get_log_prob(flat=True, discard=100))]
print(theta_max)
a, mstar, mdot = theta_max[0], theta_max[1], theta_max[2]
# a, mstar, mdot = -1, 0.15, np.log10(5e-7)
opacity_table = generate_opacity_table_x22(
    a_min=1e-6, a_max=10**(a-1), # min/max grain size
    q=-3.5, # slope for dust size distribution, dn/da ~ a^q
    dust_to_gas=0.01 # dust-to-gas ratio before sublimation
)
disk_property_table = generate_disk_property_table(opacity_table = opacity_table)
diskmodel = DiskModel(opacity_table, disk_property_table)
diskmodel.generate_disk_profile(Mstar=mstar*Msun, Mdot=10**mdot*Msun/yr, Rd=30*au, Q=1.5)
diskmodel.set_lam_obs_list([0.13])
diskmodel.generate_observed_flux(cosI=np.cos(np.deg2rad(73)))

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
DI_alma.generate_mock_observation(R=diskmodel.R, I=diskmodel.I_obs[0], cosI=np.cos(np.deg2rad(73)))
fig, ax = plt.subplots(1,3, sharex=True, sharey=False, figsize=(15,5))
fig.subplots_adjust(left=0.05, right=0.97, top=0.9, bottom=0.1, wspace=0.0, hspace=0.0)

cb68 = ax[0].imshow(DI_alma.img*1e3, cmap='jet', origin='lower', vmin=-0.1, vmax=4)
colorbar = fig.colorbar(cb68, ax=ax[0], pad=0.00, aspect=30, shrink=.98)
colorbar.set_label('Intensity (mJy/beam)')
ax[0].set_xlabel('AU', fontsize=14)
ax[0].set_yticks([0, DI_alma.img.shape[0]//2, DI_alma.img.shape[0]-1])
ax[0].set_yticklabels([-np.round(DI_alma.img_size_au), 0, np.round(DI_alma.img_size_au)])
ax[0].set_ylabel('AU', fontsize=14)
ax[0].set_title('eDisk (1.3 mm)', fontsize=14)

beam = Ellipse((130, 10), width=DI_alma.beam_min_au/DI_alma.au_per_pix, height=DI_alma.beam_maj_au/DI_alma.au_per_pix,
            angle=DI_alma.beam_pa, edgecolor='w', facecolor='w', lw=1.5, fill=True)
ax[0].add_patch(beam)


# ax[0].contour(DI_alma.img, levels=[50*40e-6]ors='black', linewidths=0.65)

model = ax[1].imshow(DI_alma.img_model*1e3, cmap='jet', origin='lower', vmin=-0.1, vmax=4)
colorbar = fig.colorbar(model, ax=ax[1], pad=0.00, aspect=30, shrink=.98)
colorbar.set_label('Intensity (mJy/beam)')
beam = Ellipse((120, 10), width=DI_alma.beam_min_au/DI_alma.au_per_pix, height=DI_alma.beam_maj_au/DI_alma.au_per_pix,
            angle=DI_alma.beam_pa, edgecolor='w', facecolor='w', lw=1.5, fill=True)
ax[1].set_xlabel('AU', fontsize=14)
ax[1].set_yticks([])
ax[1].set_title('Best-fit GIdisk model', fontsize=14)

beam = Ellipse((130, 10), width=DI_alma.beam_min_au/DI_alma.au_per_pix, height=DI_alma.beam_maj_au/DI_alma.au_per_pix,
            angle=DI_alma.beam_pa, edgecolor='w', facecolor='w', lw=1.5, fill=True)
ax[1].add_patch(beam)

residual = ax[2].imshow((DI_alma.img-DI_alma.img_model)*1e3, cmap=residual_cmp, origin='lower', vmin=-2, vmax=2)
colorbar = fig.colorbar(residual, ax=ax[2], pad=0.00, aspect=30, shrink=.98)
colorbar.set_label('Intensity (mJy/beam)')
ax[2].set_xticks([0, DI_alma.img.shape[0]//2, DI_alma.img.shape[0]-1])
ax[2].set_xticklabels([-np.round(DI_alma.img_size_au), 0, np.round(DI_alma.img_size_au)])
ax[2].set_xlabel('AU', fontsize=14)
ax[2].set_yticks([])
ax[2].set_title('Residual', fontsize=14)

beam = Ellipse((130, 10), width=DI_alma.beam_min_au/DI_alma.au_per_pix, height=DI_alma.beam_maj_au/DI_alma.au_per_pix,
            angle=DI_alma.beam_pa, edgecolor='k', facecolor='k', lw=1.5, fill=True)
ax[2].add_patch(beam)


# plt.tight_layout()
# plt.show()
plt.savefig('residual_best.pdf', transparent=True)
