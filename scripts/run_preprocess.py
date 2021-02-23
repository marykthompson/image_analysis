#run_preprocess
#run preprocessing on the acquired images, including channel split, bg subtract
#and name according to the _dapi.tif and _fish.tif conventions
#the .vsi conversion code is from: https://github.com/pskeshu/microscoper

from __future__ import unicode_literals, print_function
import numpy as np
import skimage.io as skio
import os
import os
import collections
import bioformats as bf
import javabridge as jb
import numpy as np
import tifffile as tf
import tqdm
#from .args import arguments
import xml.dom.minidom

def _init_logger():
    """This is so that Javabridge doesn't spill out a lot of DEBUG messages
    during runtime.
    From CellProfiler/python-bioformats.
    """
    rootLoggerName = jb.get_static_field("org/slf4j/Logger",
                                         "ROOT_LOGGER_NAME",
                                         "Ljava/lang/String;")

    rootLogger = jb.static_call("org/slf4j/LoggerFactory",
                                "getLogger",
                                "(Ljava/lang/String;)Lorg/slf4j/Logger;",
                                rootLoggerName)

    logLevel = jb.get_static_field("ch/qos/logback/classic/Level",
                                   "WARN",
                                   "Lch/qos/logback/classic/Level;")

    jb.call(rootLogger,
            "setLevel",
            "(Lch/qos/logback/classic/Level;)V",
            logLevel)

def get_metadata(filename):
    """Read the meta data and return the metadata object.
    """
    meta = bf.get_omexml_metadata(filename)
    metadata = bf.omexml.OMEXML(meta)
    return metadata

def get_channel(metadata, channel):
    """Return the channel name from the metadata object"""
    try:
        channel_name = metadata.image().Pixels.Channel(channel).Name
    except:
        return

    if channel_name is None:
        return
    return channel_name.replace("/", "_")

def read_image(path):
    """Reads images from the .vsi and associated files.
    Returns a dictionary with key as channel, and list
    of images as values."""
    with bf.ImageReader(path) as reader:
        # Shape of the data
        c_total = reader.rdr.getSizeC()
        z_total = reader.rdr.getSizeZ()
        t_total = reader.rdr.getSizeT()

        # Since we don't support hyperstacks yet...
        if 1 not in [z_total, t_total]:
            raise TypeError("Only 4D images are currently supported.")

        metadata = get_metadata(path)

        # This is so we can manually set a description down below.
        pbar_c = tqdm.tqdm(range(c_total))
        #store ch:np_array
        img_dict = {}
        for channel in pbar_c:
            images = []
            # Get the channel name, so we can name the file after this.
            channel_name = get_channel(metadata, channel)

            # Update the channel progress bar description with the
            # channel name.
            pbar_c.set_description(channel_name)

            for time in tqdm.tqdm(range(t_total), "T"):
                for z in tqdm.tqdm(range(z_total), "Z"):
                    image = reader.read(c=channel,
                                        z=z,
                                        t=time,
                                        rescale=False)

                    # If there's no metadata on channel name, save channels
                    # with numbers,starting from 0.
                    if channel_name is None:
                        channel_name = str(channel)

                    images.append(image)

            #save_images(np.asarray(images), channel_name, save_directory, big,
            #            save_separate)
            img_dict[channel_name] = np.asarray(images)
    #return metadata
    return metadata, img_dict

def bg_subtract(img_dict, outpath, args):
    '''
    Find all the converted .tif files in the directory, assign experiment name,
    read image, bg subtract, and save as outdir/expname_dapi.tif; outdir/expname_fish.tif
    '''
    fish_img = img_dict[str(args.fish_ch)]
    dapi_img = img_dict[str(args.dapi_ch)]
    camera_bg = args.camera_bg

    np.putmask(fish_img, fish_img < camera_bg, camera_bg)
    np.putmask(dapi_img, dapi_img < camera_bg, camera_bg)
    fish_sub = fish_img - camera_bg
    dapi_sub = dapi_img - camera_bg

    if (fish_sub < 0).any():
        print('subtraction gives vals below 0!')
    if (dapi_sub < 0).any():
        print('subtraction gives vals below 0!')

    #output as tif files. To be used for downstream segmentation.
    outfish = f'{outpath}_fish.tif'
    outdapi = f'{outpath}_dapi.tif'

    skio.imsave(outfish, fish_sub, plugin = 'tifffile')
    skio.imsave(outdapi, dapi_sub, plugin = 'tifffile')

def run(img_file, args):
    outname = os.path.basename(img_file).rstrip('.vsi')
    outpath = os.path.join(args.outdir, 'images', outname)
    metadata, img_dict = read_image(img_file)
    bg_subtract(img_dict, outpath, args)
