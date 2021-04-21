#Manage running of the bigfish pipeline
import os
import sys
import argparse
import yaml
import glob
import pandas as pd
import run_cp
import postprocess
from run_spot_quant import SpotDetector
import run_preprocess
import bioformats as bf
import javabridge as jb

def preprocess(args):
    #1) Preprocess images:
    #1a) convert from vsi to tiff
    #1b) Read tiffs and bg subtract. Put images in images/ directory
    jb.start_vm(class_path=bf.JARS, max_heap_size="2G")
    logger = run_preprocess._init_logger()
    os.makedirs(os.path.join(args.outdir, 'images'), exist_ok = True)
    vsi_files = [f.path for f in os.scandir(args.acq_dir) if f.name.endswith('.vsi')]
    #check if 'included_images or excluded images' is in args
    if hasattr(args, 'included_images'):
        imagelist = [os.path.join(args.acq_dir, i) for i in args.included_images]
        vsi_files = [i for i in vsi_files if i in imagelist]
    elif hasattr(args, 'excluded_images'):
        excluded_images = [os.path.join(args.acq_dir, i) for i in args.excluded_images]
        vsi_files = [i for i in vsi_files if i not in excluded_images]
    for i in vsi_files:
        run_preprocess.run(i, args)
    jb.kill_vm()

def run_cellpose(indir, args):
    #run cellPose segmentation on dapi images
    run_cp.run(indir, args.size_scale, args.nuc_diameter)

def bigfish(indir, args):
    #run bigfish quantification on each fov:
    infiles = os.listdir(indir)
    fov_names = [i.split('fish.tif')[0] for i in infiles if i.endswith('fish.tif')]
    fov_dict = {e:{'rna':f'{indir}/{e}fish.tif', 'nuc':f'{indir}/{e}dapi.tif'} for e in fov_names}
    for i in fov_dict:
        nuc_mask_file = os.path.join(args.outdir,
        'segmentation-results', f'{i}dapi__mask__nuclei.png')
        detector = SpotDetector(i, fov_dict[i]['rna'], fov_dict[i]['nuc'], args,
                                nuc_mask=nuc_mask_file)
        detector.run()
        detector.write_log(args)

    postprocess.summarize_experiments(args.outdir, os.path.join('results', f'{args.experiment_name}'))

def main(arglist):
    parser = argparse.ArgumentParser()
    parser.add_argument('acq_dir', help = 'folder containing acquired images')
    parser.add_argument('outdir', help = 'folder name for output directory')
    parser.add_argument('--config', help = '.yml file containing parameters for analysis')
    parser.add_argument('--preprocess', action = 'store_true', default = False, help = 'only run preprocessing')
    #If --bigfish, acq_dir should be the directory containing the prepared _fish.tif and _dapi.tif files
    parser.add_argument('--cellpose', action = 'store_true', default = False, help = 'only run cellpose')
    parser.add_argument('--bigfish', action = 'store_true', default = False, help = 'only run bigfish pipeline.')
    #parser.add_argument('--imagelist', help = 'text file containing names of images to process instead of all images')
    args = parser.parse_args()

    ##TODO##
    # - add option for no cellPose segmentation, just report foci and RNA positions
    #  independently of cell and nuclear overlap

    #get params from the config file
    with open(args.config) as file:
        config = yaml.load(file, Loader=yaml.FullLoader)

    #update args with values from the config file
    for k in config.keys():
        setattr(args, k, config[k])

    if args.preprocess:
        preprocess(args)
    elif args.cellpose:
        #run_cellpose(os.path.join(args.acq_dir, 'images'), args)
        run_cellpose(os.path.join(args.outdir, 'images'), args)
    elif args.bigfish:
        #bigfish(os.path.join(args.acq_dir, 'images'), args)
        bigfish(os.path.join(args.outdir, 'images'), args)
    else:
        preprocess(args)
        run_cellpose(os.path.join(args.outdir, 'images'), args)
        #indirectory is now the directory with preprocessed images
        bigfish(os.path.join(args.outdir, 'images'), args)

if __name__ == '__main__':
    main(sys.argv[1:])
