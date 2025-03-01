import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.modeling.models import Gaussian2D
from astropy.modeling.fitting import LevMarLSQFitter
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt
import warnings


def find_center(fname,
                x_lim=None,
                y_lim=None):
    
    if fname is None:
        print('No file to find source center!')
        return
    
    hdul = fits.open(fname)
    data = hdul[0].data
    header = hdul[0].header
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        wcs = WCS(header)
    hdul.close()
    
    # Estimate the source location (rough cropping for fitting)
    if x_lim is None:
        x_lim = [data.shape[1]//2-data.shape[1]//20, data.shape[1]//2+data.shape[1]//20]
    if y_lim is None:
        y_lim = [data.shape[0]//2-data.shape[0]//20, data.shape[0]//2+data.shape[0]//20]

            
    x_min, x_max = x_lim[0], x_lim[1]  # Adjust these based on the source
    y_min, y_max = y_lim[0], y_lim[1]
    cropped_data = data[y_min:y_max, x_min:x_max]
    
    # Create grid for fitting
    y, x = np.mgrid[:cropped_data.shape[0], :cropped_data.shape[1]]

    # Initialize a 2D Gaussian model with an initial guess
    gaussian_init = Gaussian2D(amplitude=np.max(cropped_data),
                            x_mean=(x_max - x_min) // 2,
                            y_mean=(y_max - y_min) // 2,
                            x_stddev=5, y_stddev=5)

    # Fit the Gaussian model
    fitter = LevMarLSQFitter()
    fitted_gaussian = fitter(gaussian_init, x, y, cropped_data)
    
    peak_intensity = fitted_gaussian.amplitude.value
    x_center_fit = fitted_gaussian.x_mean.value + x_min
    y_center_fit = fitted_gaussian.y_mean.value + y_min
    position_angle = np.degrees(fitted_gaussian.theta.value)
    # Convert pixel coordinates to ICRS
    center_skycoord = wcs.pixel_to_world(x_center_fit, y_center_fit)
    ra_deg, dec_deg = center_skycoord.ra.degree, center_skycoord.dec.degree
    
    return ra_deg, dec_deg, position_angle



# # Load the FITS file
# filename = '/run/media/hyp0515/fd14f880-ba6f-450f-b82d-98ba3710dc5f/backup/CB68_eDisk/CB68_SBLB_continuum_robust_0.0.image.tt0.fits'
# hdul = fits.open(filename)
# data = hdul[0].data
# header = hdul[0].header
# wcs = WCS(header)
# hdul.close()

# # Estimate the source location (rough cropping for fitting)
# x_min, x_max = 2750, 3250  # Adjust these based on the source
# y_min, y_max = 2750, 3250
# cropped_data = data[y_min:y_max, x_min:x_max]

# # Create grid for fitting
# y, x = np.mgrid[:cropped_data.shape[0], :cropped_data.shape[1]]

# # Initialize a 2D Gaussian model with an initial guess
# gaussian_init = Gaussian2D(amplitude=np.max(cropped_data),
#                            x_mean=(x_max - x_min) // 2,
#                            y_mean=(y_max - y_min) // 2,
#                            x_stddev=5, y_stddev=5)

# # Fit the Gaussian model
# fitter = LevMarLSQFitter()
# fitted_gaussian = fitter(gaussian_init, x, y, cropped_data)

# # Extract the fitted parameters
# peak_intensity = fitted_gaussian.amplitude.value
# x_center_fit = fitted_gaussian.x_mean.value + x_min
# y_center_fit = fitted_gaussian.y_mean.value + y_min
# position_angle = np.degrees(fitted_gaussian.theta.value)
# # Convert pixel coordinates to ICRS
# center_skycoord = wcs.pixel_to_world(x_center_fit, y_center_fit)

# # Print the results
# print(f"Peak Intensity: {peak_intensity:.3e} Jy/beam")
# print(f"Center Position (ICRS): {center_skycoord.to_string('hmsdms')}")
# print(f"Position Angle: {position_angle:.2f} degrees")
# # Plot the data and the fitted model
# plt.figure(figsize=(10, 8))
# plt.subplot(1, 2, 1)
# plt.imshow(cropped_data, origin='lower', cmap='viridis')
# plt.title("Cropped Data")
# plt.colorbar(label='Intensity')

# plt.subplot(1, 2, 2)
# plt.imshow(fitted_gaussian(x, y), origin='lower', cmap='viridis')
# plt.title("Fitted Gaussian")
# plt.colorbar(label='Intensity')
# plt.show()
