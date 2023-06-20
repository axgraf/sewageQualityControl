import itertools
import argparse
import numpy as np

from lib.arcgis import *
from lib.flags import Flag
from lib.utils import *
import logging

logging.basicConfig(format='%(asctime)s - %(levelname)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S',
                    stream=sys.stdout, level=logging.INFO)


class SewageQuality:

    def __init__(self, min_biomarker_threshold, minimal_number_measurements_for_mean_biomarker,
                 max_number_measurements_for_mean_biomarker):
        self.sewage_samples = None
        self.min_biomarker_threshold = min_biomarker_threshold
        self.minimal_number_measurements_for_mean_biomarker = minimal_number_measurements_for_mean_biomarker
        self.max_number_measurements_for_mean_biomarker = max_number_measurements_for_mean_biomarker
        self.biomarker_columns = ['biomarker_N1', 'biomarker_N2', 'biomarker_N3', 'biomarker_E', 'biomarker_ORF',
                                  'biomarker_RDRP']

    def __check_comments(self, sample_location, measurements_df):
        """ Flags a sewage sample if any not was given in column 'bem_lab' or 'bem_pn' """
        measurements_df[Flag.COMMENT_NOT_EMPTY.name] = \
            np.where(((measurements_df['bem_lab'].notnull()) & (measurements_df['bem_lab'] != "")) |
                     ((measurements_df['bem_pn'].notnull()) & (measurements_df['bem_pn'] != "")), True, False)
        num_flagged_rows = len(measurements_df[measurements_df[Flag.COMMENT_NOT_EMPTY.name] == True])
        logging.info("[Check comments] - [Sample location: '{}'] - {}/{} "
                     "samples were flagged due to non empty comments.".
                     format(sample_location, num_flagged_rows, measurements_df.shape[0]))

    def __biomarker_below_threshold(self, sample_location, measurements_df):
        """
        Mark biomarker values that are below given threshold or empty.
        """
        biomarkers_below_dict = dict()
        for biomarker in self.biomarker_columns:
            measurements_df[biomarker + "_below_threshold"] = np.where(
                ((measurements_df[biomarker].isnull()) | (measurements_df[biomarker] < self.min_biomarker_threshold)), True, False)
            biomarkers_below_dict[biomarker] = len(measurements_df[measurements_df[biomarker + "_below_threshold"] == True])
        logging.info("[Biomarker below threshold] - [Sample location: '{}'] - \n"
                     "\tNumber of excluded biomarkers which are below threshold or empty:\n\t\t{}\n"
                     "\tout of {} samples.".
                     format(sample_location, biomarkers_below_dict, measurements_df.shape[0]))

    def __calculate_biomarker_ratios(self, sample_location, measurements_df):
        """
        Calculate pairwise biomarker ratios only if both are above the minimal threshold level.
        In case one of the biomarker values is below the threshold the ratio is not calculated.
        """
        for biomarker1, biomarker2 in itertools.combinations(self.biomarker_columns, 2):
            measurements_df[biomarker1 + "/" + biomarker2] = \
                np.where(measurements_df[biomarker1 + "_below_threshold"] | measurements_df[biomarker2 + "_below_threshold"],
                    None,  # if one biomarker is below threshold
                    measurements_df[biomarker1] / measurements_df[biomarker2])

    def run_quality_control(self):
        """
        Main method to run the quality checks and normalization
        """
        # arcgis = Arcgis(Config(self.config))
        # self.sewage_samples = arcgis.obtain_sewage_samples()
        # Todo: switch to real data import
        import pickle
        with open('sewageData.dat', 'rb') as f:
            self.sewage_samples = pickle.load(f)

        for sample_location, measurements in self.sewage_samples.items():
            measurements = convert_sample_list2pandas(measurements)
            measurements.sort_values(by='collectionDate', ascending=False, inplace=True)  # Sort by collection date. Newest first.
            # -----------------  BIOMARKER QC -----------------------
            # 1. check for comments. Flag samples that contain any commentary
            self.__check_comments(sample_location, measurements)
            # 2. Mark biomarker values below threshold which are excluded from the analysis.
            self.__biomarker_below_threshold(sample_location, measurements)
            # 3. Calculate pairwise biomarker values if biomarkers were not marked to be below threshold.
            self.__calculate_biomarker_ratios(sample_location, measurements)
            # 4. Obtain mean biomarker ratios of biomarker pairs from last N samples
            # Todo: @Alex
            # --------------------  SUROGATVIRUS QC -------------------
            # Todo: @Lisa


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Sewage qPCR quality control",
        epilog="author: Dr. Alexander Graf (graf@genzentrum.lmu.de)")
    parser.add_argument('-b', '--biomarker_min_threshold', metavar="FLOAT", default=8, type=float,
                        help="Number of last measurements to use for calculating the mean biomarker values",
                        required=False)
    parser.add_argument('-m', '--minimal_number_measurements_for_mean_biomarker', metavar="INT", default=8, type=int,
                        help="Minimal number of measurements required for calculating mean biomarker values, otherwise this step is skipped",
                        required=False)
    parser.add_argument('-x', '--max_number_measurements_for_mean_biomarker', metavar="INT", default=50, type=int,
                        help="Number of last measurements to use for calculating the mean biomarker values",
                        required=False)
    args = parser.parse_args()
    sewageQuality = SewageQuality(args.biomarker_min_threshold,
                                  args.minimal_number_measurements_for_mean_biomarker,
                                  args.max_number_measurements_for_mean_biomarker)
    sewageQuality.run_quality_control()
