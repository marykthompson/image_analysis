'''
Create an outline file for FISHQuant starting with an ImageJ ROI file.
'''

import sys
import argparse
from datetime import datetime
from roifile import ImagejRoi

class Outline(object):
    def __init__(self, roifile, xy, z, RI, Ex, Em, NA, mscope):
        self.roifile = roifile
        self.imgname = roifile.split('.')[0]
        self.date = datetime.today().strftime('%d-%b-%Y')
        self.pix_xy = xy
        self.pix_z = z
        self.RI = RI
        self.Ex = Ex
        self.Em = Em
        self.NA = NA
        self.mscope = mscope

    def write_header(self):
        header_string = (
        f'FISH-QUANT\t\n'
        f'File-version\t3D_v1\n'
        f'RESULTS OF SPOT DETECTION PERFORMED ON {self.date} \n'
        f'COMMENT\tAutomated outline definition (batch or quick-save)\n'
        f'IMG_Raw\t{self.imgname}.tif\n'
        f'IMG_Filtered\t\n'
        f'IMG_DAPI\t\n'
        f'IMG_TS_label\t\n'
        f'FILE_settings\t\n'
        f'PARAMETERS\n'
        f'Pix-XY\tPix-Z\tRI\tEx\tEm\tNA\tType\n'
        f'{self.pix_xy}\t{self.pix_z}\t{self.RI}\t{self.Ex}\t{self.Em}\t{self.NA}\t{self.mscope}\n'
        )
        self.outfile.write(header_string)

    def write_rois(self):
        rois = ImagejRoi.fromfile(self.roifile)
        for i, r in enumerate(rois):
            self.outfile.write(f'CELL_START\tCell_{i+1}\n')
            xcoords, ycoords = r.coordinates().transpose()
            xstring = '\t'.join([str(i) for i in xcoords])
            ystring = '\t'.join([str(i) for i in ycoords])
            self.outfile.write(f'X_POS\t{xstring}\t\n')
            self.outfile.write(f'Y_POS\t{ystring}\t\n')
            self.outfile.write(f'Z_POS\t\n')
            self.outfile.write(f'CELL_END\n')

    def write(self):
        self.outfile = open(f'{self.imgname}__outline.txt', 'w')
        self.write_header()
        self.write_rois()
        self.outfile.close()

def main(arglist):
    parser = argparse.ArgumentParser()
    parser.add_argument('roifile', help = 'ImagejRoi file')
    parser.add_argument('-xy', help = 'pixels in xy')
    parser.add_argument('-z', help = 'pixels in z')
    parser.add_argument('-RI', help = 'refractive index')
    parser.add_argument('-NA', help = 'numerical aperature')
    parser.add_argument('-Em', help = 'emmission')
    parser.add_argument('-Ex', help = 'excitation')
    parser.add_argument('-microscope', default = 'confocal')
    args = parser.parse_args()

    o = Outline(args.roifile, args.xy, args.z, args.RI, args.Ex, args.Em, args.NA, args.microscope)
    o.write()

if __name__ == '__main__':
    main(sys.argv[1:])
