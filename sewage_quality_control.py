import argparse
import numpy as np
from lib.arcgis import *
from lib.flags import Flag
from lib.biomarkerQC import BiomarkerQC
import lib.utils as utils


class SewageQuality:

    def __init__(self, verbosity, quiet, min_biomarker_threshold, minimal_number_measurements_for_outlier_detection,
                 max_number_measurements_for_outlier_detection, min_number_of_valid_biomarkers):
        self.sewage_samples = None
        self.verbosity = verbosity
        self.quiet = quiet
        self.min_biomarker_threshold = min_biomarker_threshold
        self.minimal_number_measurements_for_outlier_detection = minimal_number_measurements_for_outlier_detection
        self.max_number_measurements_for_outlier_detection = max_number_measurements_for_outlier_detection
        self.min_number_of_valid_biomarkers = min_number_of_valid_biomarkers
        self.logger = utils.SewageLogger(verbosity, quiet)

    def __check_comments(self, sample_location, measurements_df):
        """ Flags a sewage sample if any not was given in column 'bem_lab' or 'bem_pn' """
        measurements_df[Flag.COMMENT_NOT_EMPTY.name] = \
            np.where(((measurements_df['bem_lab'].notnull()) & (measurements_df['bem_lab'] != "")) |
                     ((measurements_df['bem_pn'].notnull()) & (measurements_df['bem_pn'] != "")), True, False)
        num_flagged_rows = len(measurements_df[measurements_df[Flag.COMMENT_NOT_EMPTY.name] == True])
        self.logger.log.info("[Check comments] - [Sample location: '{}'] - {}/{} "
                         "samples were flagged due to non empty comments.".
                         format(sample_location, num_flagged_rows, measurements_df.shape[0]))

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
        biomarkerQC = BiomarkerQC(self.min_biomarker_threshold, self.minimal_number_measurements_for_outlier_detection,
                                  self.max_number_measurements_for_outlier_detection, self.min_number_of_valid_biomarkers)
        for sample_location, measurements in self.sewage_samples.items():
            self.logger.log.info("\n####################################################\n"
                             "\tSewage location: {} \n"
                             "####################################################".format(sample_location))
            measurements = utils.convert_sample_list2pandas(measurements)
            self.logger.log.info("Total number of measurements:\t{}".format(measurements.shape[0]))
            measurements.sort_values(by='collectionDate', ascending=True, inplace=True,
                                     ignore_index=True)  # Sort by collection date. Newest first.
            # -----------------  BIOMARKER QC -----------------------
            # 1. check for comments. Flag samples that contain any commentary.
            self.__check_comments(sample_location, measurements)
            # 2. Mark biomarker values below threshold which are excluded from the analysis.
            biomarkerQC.biomarker_below_threshold(sample_location, measurements)
            # 3. Calculate pairwise biomarker values if biomarkers were not marked to be below threshold.
            biomarkerQC.calculate_biomarker_ratios(sample_location, measurements)
            # 4. Detect outliers
            biomarkerQC.detect_outliers(sample_location, measurements)
            # 5. Select measurements with enough biomarker values which are not outliers
            biomarkerQC.select_minimal_required_biomarkers(sample_location, measurements)
            # 6.
            biomarkerQC.report_last_biomarkers_invalid(sample_location, measurements)

            # Todo: @Alex
            # --------------------  SUROGATVIRUS QC -------------------
            # Todo: @Lisa


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Sewage qPCR quality control",
        epilog="author: Dr. Alexander Graf (graf@genzentrum.lmu.de)")
    parser.add_argument('-b', '--biomarker_min_threshold', metavar="FLOAT", default=8, type=float,
                        help="Minimal biomarker threshold",
                        required=False)
    parser.add_argument('-m', '--min_number_measurements_for_outlier_detection', metavar="INT", default=9, type=int,
                        help="Minimal number of previous measurements required for outlier detection, otherwise this step is skipped",
                        required=False)
    parser.add_argument('-x', '--max_number_measurements_for_outlier_detection', metavar="INT", default=50, type=int,
                        help="Maximal number of previous measurements to use for outlier detection",
                        required=False)
    parser.add_argument('-t', '--min_number_of_valid_biomarkers', metavar="INT", default=2, type=int,
                        help="The minimal number of biomarkers required",
                        required=False)
    parser.add_argument('-v', '--verbosity', action="count", help="Increase output verbosity. Can be used multiple times.")
    parser.add_argument('-q', '--quiet', action='store_true', help="Print litte output.")
    args = parser.parse_args()
    sewageQuality = SewageQuality(args.verbosity, args.quiet, args.biomarker_min_threshold,
                                  args.min_number_measurements_for_outlier_detection,
                                  args.max_number_measurements_for_outlier_detection,
                                  args.min_number_of_valid_biomarkers)
    sewageQuality.run_quality_control()
