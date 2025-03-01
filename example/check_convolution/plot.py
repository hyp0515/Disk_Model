import numpy as np
import matplotlib.pyplot as plt
from radmc3dPy import image

conti = image.readImage(fname='./test/outfile/conti_test_scat.out')

sizeau = 1000
npix = 500
pixel_area = (sizeau/npix/140)**2
beam_axis = [0.68, 0.58]
beam_area = beam_axis[0]*beam_axis[1]*np.pi/(4*np.log(2))


conti_conv = conti.imConv(dpc=140, fwhm=beam_axis, pa=0)

fig, ax = plt.subplots(1, 2)

ax[0].imshow(conti.imageJyppix[:,:,0].T, origin='lower', cmap='hot')
ax[0].text(0.1, 0.1, f'{np.sum(conti.imageJyppix[:,:,0]):.2f}', fontsize=12, color='white')
ax[1].imshow(conti_conv.imageJyppix[:,:,0].T, origin='lower', cmap='hot')
ax[1].text(0.1, 0.1, f'{np.sum(conti_conv.imageJyppix[:,:,0])/(140**2):.4f}', fontsize=12, color='white')

print(np.max(conti_conv.imageJyppix[:,:,0])/(140**2))

plt.show()