#Post processing functions for FISH analysis
import glob
import os
import pandas as pd
import bigfish.stack as stack
import numpy as np
from skimage.measure import regionprops
from skimage.segmentation import clear_border

def avg_area_whole_objects(infiles):
    '''
    Report average area of whole objects in mask files,
    e.g. way to get average area across all nuclei in a dataset
    '''
    areas = []
    for i in infiles:
        label = stack.read_image(i)
        whole_label = clear_border(label)
        props = regionprops(whole_label)
        areas.extend([i.area for i in props])
    return np.mean(areas)

def count_nuclei(mask_file, mode='whole', avg_area=None):
    '''
    Return the nuclei count of the image according to the following modes
    whole = only count whole nuclei
    all = count all nuclei
    fractional = estimate fractional nuclei of ones touching the border
    avg_area = average area (calculated elsewhere, e.g. from a whole dataset)
    to use for calculating the fractional area
    '''
    nuc_label = stack.read_image(mask_file)
    #separate into whole nuclei and fractional nuclei
    whole_label = clear_border(nuc_label)
    nuc_props = regionprops(nuc_label)
    whole_props = regionprops(whole_label)

    if mode == 'whole':
        return len(whole_props)
    elif mode == 'all':
        return len(nuc_props)
    elif mode == 'fractional':
        whole_labels = [i.label for i in whole_props]
        if not avg_area:
            avg_area = np.mean([i.area for i in whole_props])
        touching_border = [i for i in nuc_props if i.label not in whole_labels]
        border_area = np.sum([i.area for i in touching_border])
        extra_nuclei = border_area/avg_area
        return len(whole_props) + extra_nuclei

def summarize_experiments(indir, outname, spots_ext='_spots_and_foci.npy', foci_ext='_foci.npy',
                          nuc_ext='_dapi__mask__nuclei.png', nuc_count_mode='fractional',
                          res_subdir='results', nuc_subdir='segmentation-results'):
    '''
    Summarize results of all experiments in the dataset.
    indir = path to main output folder
    outname = path to output csv file, will be within the indir
    spots_ext = extension of spots npy file (after foci assignment)
    foci_ext = extension of foci file
    nuc_ext = extension of the nuclear mask files
    res_subdir = subdirectory for the npy files
    nuc_subdir = subdirectory for the nuclear mask files
    nuc_count_mode = how to count nuclei (see count_nuclei())
    '''

    #get average nuclear area across the dataset
    search_path = f'{indir}/{nuc_subdir}/*{nuc_ext}'
    mask_files = glob.glob(f'{indir}/{nuc_subdir}/*{nuc_ext}')
    avg_nuc_area = avg_area_whole_objects(mask_files)

    #regenerate the summary df using all nuclei as the nuclear masks
    array_files = glob.glob(f'{indir}/{res_subdir}/*{spots_ext}')
    res_dict = {}
    for i in array_files:
        expname = os.path.basename(i).split(spots_ext)[0]
        res_dict[expname] = {}
        foci_file = os.path.join(indir, res_subdir, f'{expname}{foci_ext}')
        mask_file = os.path.join(indir, nuc_subdir, f'{expname}{nuc_ext}')
        nuc_count = count_nuclei(mask_file, mode='fractional', avg_area=avg_nuc_area)
        nuc_label = stack.read_image(mask_file)
        rnas = np.load(i)
        foci = np.load(foci_file)
        ts, non_ts = stack.identify_objects_in_region(nuc_label, foci, 3)
        nuc_rnas, cyt_rnas = stack.identify_objects_in_region(nuc_label, rnas, 3)

        total = len(rnas)
        ts = ts[:,3].sum()
        nonts = total - ts
        nuc = len(nuc_rnas)
        cyt = len(cyt_rnas)

        res_dict[expname]['total'] = total
        res_dict[expname]['ts'] = ts
        res_dict[expname]['non_ts'] = nonts
        res_dict[expname]['nuclear'] = nuc
        res_dict[expname]['cytoplasmic'] = cyt
        res_dict[expname]['nuc_count'] = nuc_count

    res_df = pd.DataFrame.from_dict(res_dict, orient='index')
    res_df['ts/nuc'] = res_df['ts']/res_df['nuc_count']
    res_df['non_ts/nuc'] = res_df['non_ts']/res_df['nuc_count']
    res_df.to_csv(os.path.join(indir, f'{outname}_summary.csv'))
