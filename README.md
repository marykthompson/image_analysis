# FISH quantification with BigFISH wrapper:

This workflow will analyze .vsi files in a given directory with BigFISH to
count primary and mature transcripts and nuclei.
Currently set up to analyze NMJ images that don't need cell segmentation. In
the future I will add more options for segmentation.

## Usage

#### Step 1: create a configuration file with microscope and analysis parameters
This is a yaml file. See examples in param_examples/

#### Step 2: Run the main workflow:
    python run_bigfish.py <indir> <outdir> --config <config_file>
indir = path to images  
outdir = output path

#### Additional options:
The following flags can be used to run only the specified portion of the workflow:  
--prepocess: run only the vsi conversion step  
--cellpose: run only the bigfish nuclear segmentation  
--bigfish: run only the bigfish part of the pipeline

## Functions that can be used for postprocessing other data:
The following functions are integrated into the pipeline, but might also be
useful to run on other datasets or images that were processed elsewhere.

#### Summarize results
To summarize results of all fields of view into an csv file:
```
postprocess.summarize_experiments(<indir>, <outname>)
```
#### Plot 3D representation of spots and foci
To plot a 3D view of all identified spots and foci
```
fish_plots.plot_3d_foci(<img_array>, <spots_and_foci_array>)
```
img_array: numpy array of FISH image. Can be produced with bigfish.stack.read_image().  
spots_and_foci: numpy array of spots and assigned foci. This is produced by the
bigfish.detection.detect_foci() function.
