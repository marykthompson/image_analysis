#Run segmentation with CellPose
import segwrap
from segwrap import  utils_segmentation
from segwrap import utils_cellpose
from pathlib import Path
import os

def max_proj(img_dir, outdir):
    #1) Create the max projections as input for CellPose
    # Parameters
    path_process = Path(r'%s' % img_dir) # For example data: \example_data\acquisition
    path_save = Path(r'%s' % outdir)
    channel_ident = 'dapi'                           # Identifier of channel that should be pre-processed
    img_ext = '.tif'                                 # Extension of images that should be processed
    projection_type = 'max'                          # Projection type (mean, max, indiv)
    search_recursive = False                         # Search folder recursively for data?

    # Call pre-processing function
    utils_segmentation.folder_prepare_prediction(
                path_process=path_process,
                channel_ident=channel_ident,
                img_ext=img_ext,
                path_save=path_save,
                projection_type=projection_type,
                search_recursive=search_recursive)

def segment(img_dir, outdir, size_scale, nuc_diameter):
    #2) Run CellPose segmentation
    # Parameters
    path_scan = Path(r'%s' % img_dir)            # For example data: example_data\analysis\segmentation-input
    path_save = Path(r'%s' % outdir)
    obj_name='nuclei',                                 # Name of object that should be segmented
    str_channel = 'dapi'                                # Identifier of channel for nuclear segmentation
    img_ext = '.png'                                   # Extension of images to be segmented
    new_size = (size_scale,)                                   # Size of images (tuple with 1 element for scaling factor, tuple with 2 elements for new size)
    obj_size = nuc_diameter                                   # Typical size (diameter) of nuclei (in px, of the rescaled image)
    model_type = 'nuclei'                            # Can be 'nuclei' or 'cyto', for densely packed nuclei 'cyto' might work well.

    # Call segmentation function
    utils_cellpose.segment_obj_indiv(
                                    path_scan=path_scan,
                                    obj_name='nuclei',
                                    str_channel=str_channel,
                                    img_ext=img_ext,
                                    new_size=new_size,
                                    obj_size= obj_size,
                                    model_type=model_type,
                                    path_save=path_save)

def run(img_dir, size_scale, nuc_diameter):
    '''Run the CellPose segmentation'''
    outdir = img_dir.split(os.path.basename(os.path.normpath(img_dir)))[0]
    seg_indir = os.path.join(outdir, 'segmentation-input')
    seg_outdir = os.path.join(outdir, 'segmentation-results')
    #take /images/ off the end of the path
    max_proj(img_dir, seg_indir)
    segment(seg_indir, seg_outdir, size_scale, nuc_diameter)
