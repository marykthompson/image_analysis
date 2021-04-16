import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np

# Plots for visualizing bigfish analysis
def plot_3d_foci(img, spots_and_foci, nuc_mask=None, voxel_size_yx=110,
                 voxel_size_z=200, s=1, figsize=(8,8), fmt='png', dpi=300,
                 plot_img=True, vmax=None, outname='detected_foci'):
    '''
    Plot the individual spots and the called foci in 3D.
    img = the image array (needed for knowing the dimensions)
    spots_and_foci = the output of Bigfish detection.detect_foci()
    nuc_mask = the mask output by cellpose, e.g. dapi__mask__nuclei.png
    voxel_size_yx, z = voxel size, in nm, for relative scaling of the axes.
    s = size of points to plot
    fmt = figure format, i.e. png, svg
    plot_img = whether to include a max projection plot of the img in bottom right
    vmax = max intensity to scale the image color map. It is useful to put this
    approximately 2-fold higher than single RNAs, so that they will be easily visible.
    ##TODO: would be useful to plot nuclear mask (in 2D or 3D) if available
    '''
    #index positions in array
    x_i = 2
    y_i = 1
    z_i = 0
    fig = plt.figure(constrained_layout=False, figsize = figsize)

    #calculate relative voxel size for scaling
    voxel_s = voxel_size_yx/voxel_size_z
    spec = gridspec.GridSpec(ncols=2, nrows=2, figure=fig, width_ratios=[voxel_s, 1], height_ratios=[voxel_s, 1])
    ax_xy = fig.add_subplot(spec[0])
    ax_xz = fig.add_subplot(spec[2])
    ax_yz = fig.add_subplot(spec[1])
    ax_img = fig.add_subplot(spec[3])

    labels = spots_and_foci[:,-1]

    # Black removed and is used for noise instead.
    unique_labels = set(labels)
    unique_labels.remove(-1)
    colors = [plt.cm.Spectral(each)
              for each in np.linspace(0, 1, len(unique_labels))]

    #plot non-clustered spots first, so they will fall behind
    bg_mask = (labels == -1)
    d = spots_and_foci[bg_mask]
    ax_xy.scatter(d[:, x_i], d[:, y_i], color='lightgrey', s=s)
    ax_xz.scatter(d[:, x_i], d[:, z_i], color='lightgrey', s=s)
    ax_yz.scatter(d[:, z_i], d[:, y_i], color='lightgrey', s=s)
    if plot_img:
        img_mip = img.max(axis = 0)
        ax_img.imshow(img_mip, vmin=0, vmax=vmax)
        ax_img.xaxis.set_ticks([])
        ax_img.yaxis.set_ticks([])
        ax_img.set_title('original image')

    #plot each cluster with its own color
    for k, col in zip(unique_labels, colors):
        class_member_mask = (labels == k)
        d = spots_and_foci[class_member_mask]
        ax_xy.scatter(d[:, x_i], d[:, y_i], color=tuple(col), s=s)
        ax_xz.scatter(d[:, x_i], d[:, z_i], color=tuple(col), s=s)
        ax_yz.scatter(d[:, z_i], d[:, y_i], color=tuple(col), s=s)

    ax_xy.xaxis.set_ticks([])
    ax_xy.yaxis.set_ticks([])
    zmax = img.shape[z_i]
    ymax = img.shape[y_i]
    xmax = img.shape[x_i]
    ax_yz.yaxis.set_ticks_position('right')
    ax_yz.yaxis.set_label_position('right')
    ax_yz.xaxis.set_ticks_position('top')
    ax_yz.xaxis.set_label_position('top')
    ax_yz.set_xlabel('z position')
    ax_yz.set_ylabel('y position')
    ax_xz.set_ylabel('z position')
    ax_xz.set_xlabel('x position')
    ax_xy.set_title('detected spots')
    #yaxis needs to be inverted in order to display in the correct way with 0 at the top left
    ax_xy.set_ylim(ymax, 0)
    ax_xy.set_xlim(0, xmax)
    ax_yz.set_ylim(ymax, 0)
    ax_yz.set_xlim(0, zmax)
    ax_xz.set_xlim(0, xmax)
    ax_xz.set_ylim(zmax, 0)
    plt.savefig(f'{outname}.{fmt}', dpi=dpi)
