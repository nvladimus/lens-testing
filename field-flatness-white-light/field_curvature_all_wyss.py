import numpy as np
import matplotlib.pyplot as plt
import tifffile as tif
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import os
from pathlib import Path


def contrast(roi):
    """Compute the contrast value, (max-min)/(max+min), from the image roi"""
    mini = np.percentile(roi, 1)
    maxi = np.percentile(roi, 99)
    contrast = (maxi - mini) / (maxi + mini)
    return contrast


def main():
    datafolder = "Y:\\Alice\\9007_CBT_WYSS_LIGHTSHEET\\calibration_test\\mesoSPIM_detection_characterisation"
    LENS_NAME = "Olympus MVPLAPO-1x"
    datafiles = {
        "1x": datafolder + "\\mesoLens_1x.tif",
        "1.25x": datafolder + "\\mesoLens_1_25x.tif",
        "2x": datafolder + "\\mesoLens_2x.tif",
        "3.2x": datafolder + "\\mesoLens_3_2x.tif",
        "4x": datafolder + "\\mesoLens_4x.tif",
        "5x": datafolder + "\\mesoLens_5x.tif",
        "6.3x": datafolder + "\\mesoLens_6_3x.tif"
    }

    save_figs_folder = "C:\\Users\\YCABARA\\ALICe\\Quality\\Lens-test-results\\"
    f_step_um = 10
    CAMERA = "Hamamatsu Orca Flash 4.3"
    sensor_dim_mm = (13.312, 13.312)
    N_ROIs_H, N_ROIs_W = 16, 16

    zoom_keys = datafiles.keys()  # Get all the zoom keys

    for zoom_key in zoom_keys:
        MAG = float(zoom_key[:-1])  # Effective lens magnification

        # Read the image for the current zoom key
        img = tif.imread(datafiles[zoom_key])
        im_z, im_h, im_w = img.shape
        assert img.max() < 65535, "Error: image is saturated"

        roi_h, roi_w = int(im_h / N_ROIs_H), int(im_w / N_ROIs_W)

        contrast_table = np.empty((im_z, N_ROIs_H, N_ROIs_W))
        for f in range(im_z):
            for j in range(N_ROIs_H):
                for i in range(N_ROIs_W):
                    roi = img[f, j * roi_h:(j + 1) * roi_h, i * roi_w:(i + 1) * roi_w]
                    contrast_table[f, j, i] = contrast(roi)

        # Plotting Contrast Maps
        fig = plt.figure(figsize=(12, 12))
        fig.suptitle(f'Contrast maps: {LENS_NAME} zoom {zoom_key}', fontsize=20)

        ax1 = plt.subplot(2, 2, 1)
        ax1.imshow(contrast_table[z_contrast_max, :, :], vmin=cmin, vmax=cmax,
                   extent=[0, x_range[-1], 0, y_range[-1]],
                   aspect='auto', cmap='jet', interpolation='bicubic')
        ax1.set_title(f"Highest-contrast plane, F pos. {z_contrast_max * f_step_um} um", fontsize=16)
        ax1.set_xlabel("FOV_X, mm", fontsize=16)
        ax1.set_ylabel("FOV_Y, mm", fontsize=16)
        ax1.axhline(y_range[-1] / 2, c="gray", linewidth=1)
        ax1.axvline(x_range[-1] / 2, c="gray", linewidth=1)

        ax3 = plt.subplot(2, 2, 2)
        conmap3 = ax3.imshow(contrast_yz.T, vmin=cmin, vmax=cmax,
                             extent=[0, f_step_um * (im_z - 1), 0, y_range[-1]],
                             aspect='auto', cmap='jet', interpolation='bicubic')
        ax3.set_title("Contrast across FOV Y-coordinate", fontsize=16)
        ax3.set_xlabel("Focus position F, um", fontsize=16)
        ax3.set_ylabel("FOV_Y, mm", fontsize=16)
        ax3.grid(True)

        ax4 = plt.subplot(2, 2, 3)
        conmap1 = ax4.imshow(contrast_xz, vmin=cmin, vmax=cmax,
                             extent=[0, x_range[-1], 0, f_step_um * (im_z - 1)],
                             aspect='auto', cmap='jet', interpolation='bicubic')
        ax4.set_title("Contrast across FOV X-coordinate", fontsize=16)
        ax4.set_ylabel("Focus position F, um", fontsize=16)
        ax4.set_xlabel("FOV_X, mm", fontsize=16)
        ax4.grid(True)

        ax0 = plt.subplot(2, 2, 4)
        conmap0 = ax0.imshow(contrast_table.max(axis=0), vmin=cmin, vmax=cmax,
                             extent=[0, x_range[-1], 0, y_range[-1]],
                             aspect='auto', cmap='jet', interpolation='bicubic')
        ax0.set_title("Maximum contrast along F", fontsize=16)
        ax0.set_xlabel("FOV_X, mm", fontsize=16)
        ax0.set_ylabel("FOV_Y, mm", fontsize=16)

        axins = inset_axes(ax3,
                           width="5%",  # width = 5% of parent_bbox width
                           height="50%",  # height : 50%
                           loc='lower left',
                           bbox_to_anchor=(1.05, 0., 1, 1),
                           bbox_transform=ax3.transAxes,
                           borderpad=0,
                           )
        plt.colorbar(conmap0, cax=axins)

        # Save contrast maps as PNG files
        fig.savefig(save_figs_folder + zoom_key + "_contrast_map.png")

        # Save contrast table as TIFF file
        tif.imwrite(save_figs_folder + zoom_key + "_contrast.tiff", contrast_table)

        # Perform other calculations, saving, etc.
        # ...

        # Display or print results
        # ...


if __name__ == "__main__":
    main()
