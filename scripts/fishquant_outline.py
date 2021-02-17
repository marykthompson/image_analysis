'''
Create an outline file for FISHQuant starting with an ImageJ ROI file.
This version can accept either roifile or tif file, or a directory containing them.
If roifile, get the outlines from the ROIs.
If tif, take the whole image as the outline.
'''

import sys
import argparse
import os
from datetime import datetime
from roifile import ImagejRoi
from PIL import Image

class Outline(object):
    def __init__(self, input, xy, z, RI, Ex, Em, NA, mscope):
        self.input = input
        if os.path.isdir(self.input):
            self.files = [os.path.join(self.input, i) for i in os.listdir(self.input) \
             if (not i.startswith('.')) and (not i.endswith('outline.txt'))]
        else:
            self.files = [self.input]
        self.extension = self.files[0].split('.')[1]
        if self.extension == 'zip':
            self.format = 'roi'
        elif self.extension == 'tif':
            self.format = 'tif'
        else:
            print('Please provide either a zip of ROIs or a tif.')
            sys.exit()
        self.date = datetime.today().strftime('%d-%b-%Y')
        self.pix_xy = xy
        self.pix_z = z
        self.RI = RI
        self.Ex = Ex
        self.Em = Em
        self.NA = NA
        self.mscope = mscope

    def write_header(self, imgname):
        header_string = (
        f'FISH-QUANT\t\n'
        f'File-version\t3D_v1\n'
        f'RESULTS OF SPOT DETECTION PERFORMED ON {self.date} \n'
        f'COMMENT\tAutomated outline definition (batch or quick-save)\n'
        f'IMG_Raw\t{os.path.basename(imgname)}.tif\n'
        f'IMG_Filtered\t\n'
        f'IMG_DAPI\t\n'
        f'IMG_TS_label\t\n'
        f'FILE_settings\t\n'
        f'PARAMETERS\n'
        f'Pix-XY\tPix-Z\tRI\tEx\tEm\tNA\tType\n'
        f'{self.pix_xy}\t{self.pix_z}\t{self.RI}\t{self.Ex}\t{self.Em}\t{self.NA}\t{self.mscope}\n'
        )
        self.outfile.write(header_string)

    def write_rois(self, img):
        rois = ImagejRoi.fromfile(img)
        for i, r in enumerate(rois):
            self.outfile.write(f'CELL_START\tCell_{i+1}\n')
            xcoords, ycoords = r.coordinates().transpose()
            xstring = '\t'.join([str(i) for i in xcoords])
            ystring = '\t'.join([str(i) for i in ycoords])
            self.outfile.write(f'X_POS\t{xstring}\t\n')
            self.outfile.write(f'Y_POS\t{ystring}\t\n')
            self.outfile.write(f'Z_POS\t\n')
            self.outfile.write(f'CELL_END\n')

    def write_tifbox(self, img):
        width, height = Image.open(img).size
        xstring = '\t'.join(['1', '1', str(width), str(width)])
        ystring = '\t'.join(['1', str(height), str(height), '1'])
        self.outfile.write(f'CELL_START\tEntireImage\n')
        self.outfile.write(f'X_POS\t{xstring}\t\n')
        self.outfile.write(f'Y_POS\t{ystring}\t\n')
        self.outfile.write(f'Z_POS\t\n')
        self.outfile.write(f'CELL_END\n')

    def write(self):
        for img in self.files:
            imgname = img.split('.')[0]
            self.outfile = open(f'{imgname}__outline.txt', 'w')
            self.write_header(imgname)
            if self.format == 'roi':
                self.write_rois(img)
            elif self.format == 'tif':
                self.write_tifbox(img)
            self.outfile.close()

def main(arglist):
    parser = argparse.ArgumentParser()
    parser.add_argument('input', help = 'ImagejRoi file (should be .zip) or .tif \
    or can be a directory containing these file types')
    parser.add_argument('-xy', help = 'pixels in xy')
    parser.add_argument('-z', help = 'pixels in z')
    parser.add_argument('-RI', help = 'refractive index')
    parser.add_argument('-NA', help = 'numerical aperature')
    parser.add_argument('-Em', help = 'emmission')
    parser.add_argument('-Ex', help = 'excitation')
    parser.add_argument('-microscope', default = 'confocal')
    args = parser.parse_args()

    o = Outline(args.input, args.xy, args.z, args.RI, args.Ex, args.Em, args.NA, args.microscope)
    o.write()

if __name__ == '__main__':
    main(sys.argv[1:])
