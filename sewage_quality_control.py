import os
import argparse
import numpy as np
from lib.arcgis import *
import pandas as pd
from lib.flags import Flag
from lib.biomarkerQC import BiomarkerQC
import lib.utils as utils


class SewageQuality:

    def __init__(self, output_folder, verbosity, quiet, outlier_statistics, min_biomarker_threshold, minimal_number_measurements_for_outlier_detection,
                 max_number_measurements_for_outlier_detection, min_number_of_valid_biomarkers_ratios, report_outliers, periode_month_surrogatevirus):
        self.sewage_samples = None
        self.output_folder = output_folder
        self.verbosity = verbosity
        self.quiet = quiet
        self.outlier_statistics = outlier_statistics
        self.min_biomarker_threshold = min_biomarker_threshold
        self.minimal_number_measurements_for_outlier_detection = minimal_number_measurements_for_outlier_detection
        self.max_number_measurements_for_outlier_detection = max_number_measurements_for_outlier_detection
        self.min_number_of_valid_biomarkers_ratios = min_number_of_valid_biomarkers_ratios
        self.periode_month_surrogatevirus = periode_month_surrogatevirus
        self.report_outliers = report_outliers
        self.logger = utils.SewageLogger(verbosity, quiet)
        self.__setup()

    def __setup(self):
        if not os.path.exists(self.output_folder):
            os.makedirs(self.output_folder)

    def __check_comments(self, sample_location, measurements_df: pd.DataFrame):
        """ Flags a sewage sample if any text is stored in column 'bem_lab' or 'bem_pn' """
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
        #arcgis = Arcgis(Config(self.config))
        #self.sewage_samples = arcgis.obtain_sewage_samples()
        # Todo: switch to real data import
        import pickle
        with open('sewageData.dat', 'rb') as f:
            self.sewage_samples = pickle.load(f)
        with open('sewagePlantData.dat', 'rb') as f:
            self.sewage_plants2trockenwetterabfluss = pickle.load(f)
        biomarkerQC = BiomarkerQC(self.output_folder, self.outlier_statistics, self.min_biomarker_threshold,
                                  self.minimal_number_measurements_for_outlier_detection,
                                  self.max_number_measurements_for_outlier_detection,
                                  self.min_number_of_valid_biomarkers_ratios,
                                  self.report_outliers)
        for idx, (sample_location, measurements) in enumerate(self.sewage_samples.items()):
            if idx == 0:
                continue   # skip first sewage location for testing  #Todo: remove before production
            self.logger.log.info("\n####################################################\n"
                             "\tSewage location: {} \n"
                             "####################################################".format(sample_location))
            measurements = utils.convert_sample_list2pandas(measurements)
            self.logger.log.info("Total number of measurements:\t{}".format(measurements.shape[0]))
            measurements['collectionDate'] = pd.to_datetime(measurements['collectionDate'], format="%Y-%m-%d")
            measurements.sort_values(by='collectionDate', ascending=True, inplace=True,
                                     ignore_index=True)  # Sort by collection date. Newest first. Ascending important!
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
            biomarkerQC.filter_required_biomarkers(sample_location, measurements)
            # 6.
            biomarkerQC.report_last_biomarkers_invalid(sample_location, measurements)

            # Todo: @Alex
            # --------------------  SUROGATVIRUS QC -------------------
            # Todo: @Lisa


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Sewage qPCR quality control",
        usage='use "python3 %(prog)s --help" for more information',
        epilog="author: Dr. Alexander Graf (graf@genzentrum.lmu.de)", formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-o', '--output_folder', metavar="FOLDER", default="sewage_qc", type=str,
                        help="Specifiy output folder. (default folder: 'sewage_qc')",
                        required=False)
    parser.add_argument('-v', '--verbosity', action="count", help="Increase output verbosity.")
    parser.add_argument('-q', '--quiet', action='store_true', help="Print litte output.")

    biomarker_qc_group = parser.add_argument_group("Biomarker quality control")
    biomarker_qc_group.add_argument('-s', '--outlier_statistics', metavar="METHOD", default="all", nargs='*',
                        help=("Which outlier detection methods should be used? Multiple selections allowed. (default: 'all')\n"
                              "Possible choices are : [lof, rf, iqr, zscore, all]\n"
                              "To select 'rf' and 'iqr' use: --outlier_statistics rf iqr \n"
                              "\tlof = local outlier factor\n"
                              "\trf = random forest\n"
                              "\tiqr = interquartile range\n"
                              "\tzscore = modified z-score\n"
                              "\tall = use all methods\n") ,
                        choices=["lof", "rf", "iqr", "zscore", "all"],
                        required=False)
    biomarker_qc_group.add_argument('-b', '--biomarker_min_threshold', metavar="FLOAT", default=8, type=float,
                        help="Minimal biomarker threshold. (default: 8)",
                        required=False)
    biomarker_qc_group.add_argument('-m', '--min_number_measurements_for_outlier_detection', metavar="INT", default=9, type=int,
                        help="Minimal number of previous measurements required for outlier detection, otherwise this step is skipped. (default: 9)",
                        required=False)
    biomarker_qc_group.add_argument('-x', '--max_number_measurements_for_outlier_detection', metavar="INT", default=50, type=int,
                        help="Maximal number of previous measurements to use for outlier detection. (default: 50)",
                        required=False)
    biomarker_qc_group.add_argument('-t', '--min_number_of_valid_biomarkers_ratios', metavar="INT", default=1, type=int,
                        help="The minimal number of biomarkers ratios required after outlier detection. (default: 1)",
                        required=False)
    biomarker_qc_group.add_argument('-r', '--report_outliers', metavar="INT", default=2, type=int,
                        help="The number of outliers identified in the last N consecutive biomarker ratios that trigger a report. (default: 2)",
                        required=False)

    surrogatevirus_group = parser.add_argument_group("Surrogate virus quality control")
    surrogatevirus_group.add_argument('-d', '--periode_month_surrogatevirus', metavar="INT", default=4, type=int,
                        help="The periode of time (month) taken into account for surrogatevirus outliers",
                        required=False)

    args = parser.parse_args()
    if args.min_number_of_valid_biomarkers_ratios < 1:
        parser.error("The minimal number of required biomarker ratios must be >= 1 (--min_number_of_valid_biomarkers_ratios).")

    sewageQuality = SewageQuality(args.output_folder, args.verbosity, args.quiet, args.outlier_statistics, args.biomarker_min_threshold,
                                  args.min_number_measurements_for_outlier_detection,
                                  args.max_number_measurements_for_outlier_detection,
                                  args.min_number_of_valid_biomarkers_ratios,
                                  args.report_outliers, args.periode_month_surrogatevirus)
    sewageQuality.run_quality_control()

