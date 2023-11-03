# Utility functions for lens testing
# Nikita Vladimirov, 2023
# GPL-3 license

import matplotlib.pyplot as plt
import numpy as np

def wavelength_to_rgb(wavelength):
    """A crude mapping of nm to RGB space for plotting purposes"""
    # Define the wavelength range (e.g., 400-700 nm)
    min_wavelength, max_wavelength = 400, 720

    # Normalize the wavelength to the 0-1 range
    normalized_wavelength = (wavelength - min_wavelength) / (max_wavelength - min_wavelength)

    # Use a colormap (e.g., jet) to map the normalized wavelength to RGB
    colormap = plt.get_cmap('turbo') # or 'jet'
    rgb = colormap(normalized_wavelength)

    return rgb

def contrast(roi):
    """Compute the contrast value, (max-min)/(max+min), from the image roi"""
    mini = np.percentile(roi, 1)
    maxi = np.percentile(roi, 99)
    contrast = (maxi - mini) / (maxi + mini)
    return contrast

