# Run spot detection and ts quantification
import bigfish
import bigfish.stack as stack
import bigfish.segmentation as segmentation
import bigfish.plot as plot
import bigfish.detection as detection
import numpy as np
import pandas as pd
import os

def calculate_psf(voxel_size_z, voxel_size_xy, Ex, Em, NA, RI, microscope):
    '''
    Use the formula implemented in Matlab version (sigma_PSF_BoZhang_v1)
    to calculate the theoretical PSF size.
    '''
    if microscope == 'widefield':
        sigma_xy = 0.225*Em/NA
        sigma_z = 0.78*RI*Em/(NA**2)
    elif microscope in {'confocal', 'nipkow'}:
        sigma_xy = 0.225/NA*Ex*Em/np.sqrt(Ex**2 + Em**2)
        sigma_z = 0.78*RI/NA**2*Ex*Em/np.sqrt(Ex**2 + Em**2)
        #original matlab code commented below:
        #widefield:
        #sigma_xy = 0.225 * lambda_em / NA ;
        #sigma_z  = 0.78 * RI * lambda_em / (NA*NA) ;
        #confocal/nipkow
        #sigma_xy = 0.225 / NA * lambda_ex * lambda_em / sqrt( lambda_ex^2 + lambda_em^2 ) ;
        #sigma_z =  0.78 * RI / (NA^2) *  lambda_ex * lambda_em / sqrt( lambda_ex^2 + lambda_em^2 ) ;

    else:
        print(f'microscope={microscope} is not a valid option')
        sys.exit()

    return sigma_z, sigma_xy

class SpotDetector(object):

    def __init__(self, fov_name, fish_img, dapi_img, args, nuc_mask=None, cell_mask=None):
        self.voxel_size_z = args.voxel_size_z
        self.voxel_size_xy = args.voxel_size_xy
        self.nuc_mask = nuc_mask
        self.cell_mask = cell_mask
        self.psf_z, self.psf_xy = calculate_psf(self.voxel_size_z, self.voxel_size_xy,
                                                args.Ex, args.Em, args.NA, args.RI,
                                                args.microscope)
        self.rna = stack.read_image(fish_img)
        self.nuc = stack.read_image(dapi_img)
        self.rna_mip = stack.maximum_projection(self.rna)
        self.fov_name = fov_name
        self.cluster_radius = args.cluster_radius
        self.nb_in_cluster = args.nb_in_cluster
        self.plotdir = os.path.join(args.outdir, 'plots')
        self.datadir = os.path.join(args.outdir, 'results')
        os.makedirs(self.plotdir, exist_ok = True)
        os.makedirs(self.datadir, exist_ok = True)

    def detect_spots(self):
        # sigma
        sigma_z, sigma_xy, sigma_xy = detection.get_sigma(self.voxel_size_z,
                                                          self.voxel_size_xy,
                                                          self.psf_z, self.psf_xy)

        sigma = (sigma_z, sigma_xy, sigma_xy)

        # LoG filter
        rna_log = stack.log_filter(self.rna, sigma)

        # local maximum detection
        mask = detection.local_maximum_detection(rna_log, min_distance=sigma)

        # thresholding
        threshold = detection.automated_threshold_setting(rna_log, mask)
        spots, _ = detection.spots_thresholding(rna_log, mask, threshold)

        #note radius
        (radius_z, radius_xy, radius_xy) = detection.get_radius(self.voxel_size_z,
                                                                self.voxel_size_xy,
                                                                self.psf_z, self.psf_xy)
        self.sigma_z = sigma_z
        self.sigma_xy = sigma_xy
        self.radius_z = radius_z
        self.radius_xy = radius_xy
        self.threshold = threshold
        self.spots = spots
        plot.plot_detection(self.rna_mip, self.spots, radius=self.radius_xy,
                            framesize=(10, 8), contrast=True, show=False,
                            path_output=f'{self.plotdir}/{self.fov_name}spots',
                            ext='png')

    def decompose_clusters(self):
        #Cluster detection: try to decompose clusters of spots into individual foci.
        # alpha impacts the number of spots per cluster
        # beta impacts the number of detected clusters
        self.spots_post_decomposition, self.clusters, self.reference_spot = \
        detection.decompose_cluster(self.rna, self.spots, self.voxel_size_z,
                                    self.voxel_size_xy, self.psf_z, self.psf_xy,
                                    alpha=0.7, beta=1)

        plot.plot_reference_spot(self.reference_spot,
                                 path_output=f'{self.plotdir}/{self.fov_name}refspot',
                                 ext='png', contrast=True, show=False)

        plot.plot_detection(self.rna_mip, self.spots_post_decomposition,
                            radius=self.radius_xy, framesize=(10, 8),
                            contrast=True, show=False,
                            path_output=f'{self.plotdir}/{self.fov_name}clusters',
                            ext='png')

    def detect_foci(self):
        #Clusters that are big enough can be considered foci
        self.spots_post_clustering, self.foci = detection.detect_foci(
                                                self.spots_post_decomposition,
                                                self.voxel_size_z, self.voxel_size_xy,
                                                self.cluster_radius, self.nb_in_cluster)

        plot.plot_detection(self.rna_mip, [self.spots_post_decomposition, self.foci[:, :3]],
                            path_output = f'{self.plotdir}/{self.fov_name}foci',
                            ext='png', shape=["circle", "polygon"],
                            radius=[self.radius_xy, self.radius_xy*2],
                            color=["red", "blue"], linewidth=[1, 2],
                            fill=[False, True], framesize=(10, 8), contrast=True,
                            show=False)

    def assign_foci_ts(self):
        '''
        Split foci into foci that overlap the nucleus (ts)
        and foci which do not overlap.
        Extract the results for the FOV.
        '''
        nuc_label = stack.read_image(self.nuc_mask)
        if self.cell_mask is None:
            cell_label = np.ones(nuc_label.shape, dtype = 'int64')

        self.spots_no_ts, self.non_ts_foci, self.ts = stack.remove_transcription_site(
                                                      self.spots_post_clustering,
                                                      self.foci, nuc_label, ndim=3)

        image_contrasted = stack.rescale(self.rna, channel_to_stretch=0)
        image_contrasted = stack.maximum_projection(image_contrasted)
        self.nuc_mip = stack.maximum_projection(self.nuc)

        #Get results for field of view
        self.fov_results = stack.extract_cell(cell_label=cell_label, ndim=3,
                                              nuc_label=nuc_label,
                                              rna_coord=self.spots_no_ts,
                                              others_coord={"foci": self.foci,
                                              "transcription_site": self.ts},
                                              image=image_contrasted,
                                              others_image={"dapi": self.nuc_mip,
                                              "smfish": self.rna_mip},
                                              remove_cropped_cell=False,
                                              check_nuc_in_cell=False)

        print("number of cells identified: {0}".format(len(self.fov_results)))

    def plot_fov(self):
        '''
        Plot detected cells, spots, nuclei in the FOV.
        '''
        for i, cell_results in enumerate(self.fov_results):
            print("cell {0}".format(i))

            # get cell results
            cell_mask = cell_results["cell_mask"]
            cell_coord = cell_results["cell_coord"]
            nuc_mask = cell_results["nuc_mask"]
            nuc_coord = cell_results["nuc_coord"]
            rna_coord = cell_results["rna_coord"]
            foci_coord = cell_results["foci"]
            ts_coord = cell_results["transcription_site"]
            image_contrasted = cell_results["image"]
            print("\r number of rna {0}".format(len(rna_coord)))
            print("\r number of foci {0}".format(len(foci_coord)))
            print("\r number of transcription sites {0}".format(len(ts_coord)))

            # It only plots the first nucleus for multinucleate cells.
            plot.plot_cell(ndim=3, cell_coord=cell_coord, nuc_coord=nuc_coord,
                          rna_coord=rna_coord, foci_coord=foci_coord,
                          other_coord=ts_coord, image=image_contrasted,
                          cell_mask=cell_mask, nuc_mask=nuc_mask,
                          title="Cell {0}".format(i), framesize=(12, 10),
                          path_output=f'{self.plotdir}/{self.fov_name}cell_{i}',
                          ext='png', show=False)

    def summarize_fov(self):
        df = stack.summarize_extraction_results(self.fov_results, ndim=3)
        #add the missing nb_nascent RNAs to the df:
        for i, cell_results in enumerate(self.fov_results):
            cell_id = cell_results['cell_id']
            df.loc[df['cell_id'] == cell_id, 'nb_nascent'] = cell_results['transcription_site'][:, 3].sum()
        df.to_csv(f'{self.datadir}/{self.fov_name}results.csv')

    def save_arrays(self):
        '''
        Save arrays:
        To save for primary analysis: spots_no_ts, ts
        To save if using different segmentation: spots_post_clustering, foci
        '''
        stack.save_array(self.spots_no_ts, f'{self.datadir}/{self.fov_name}spots_no_ts.npy')
        stack.save_array(self.ts, f'{self.datadir}/{self.fov_name}ts.npy')
        stack.save_array(self.spots_post_clustering, f'{self.datadir}/{self.fov_name}spots_post_clustering.npy')
        stack.save_array(self.foci, f'{self.datadir}/{self.fov_name}foci.npy')

    def write_log(self, args):
        '''
        Write log file containing values from config, as well as calculated
        values.
        '''
        log_dict = vars(args)
        atts = ['threshold', 'radius_z', 'radius_xy', 'psf_z', 'psf_xy',
                'sigma_z', 'sigma_xy']
        for i in atts:
            log_dict[i] = getattr(self, i)

        order = ['dapi_ch', 'fish_ch', 'camera_bg', 'size_scale', 'nuc_diameter',
        'voxel_size_xy', 'voxel_size_z', 'NA', 'RI', 'Ex', 'Em', 'microscope',
        'cluster_radius', 'nb_in_cluster']
        order.extend(atts)

        df = pd.DataFrame.from_dict(log_dict, orient='index', columns=['value'])
        df.index.name = 'parameter'
        df.reindex(order).to_csv(f'{self.datadir}/{self.fov_name}log.csv')

    def run(self):
        self.detect_spots()
        self.detect_spots()
        self.decompose_clusters()
        self.detect_foci()
        self.assign_foci_ts()
        self.plot_fov()
        self.summarize_fov()
        self.save_arrays()
