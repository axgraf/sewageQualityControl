import os
import itertools
import argparse
import numpy as np
from lib.arcgis import *
from lib.constant import *
from lib.biomarkerQC import BiomarkerQC
from lib.water_quality import WaterQuality
import lib.utils as utils


class SewageQuality:

    def __init__(self, output_folder, verbosity, quiet, biomarker_outlier_statistics, min_biomarker_threshold,
                 min_number_biomarkers_for_outlier_detection,
                 max_number_biomarkers_for_outlier_detection, report_number_of_biomarker_outlier, periode_month_surrogatevirus,
                 water_quality_number_of_last_month, min_number_of_last_measurements_for_water_qc, water_qc_outlier_statistics):
        self.sewage_samples = None
        self.output_folder = output_folder
        self.verbosity = verbosity
        self.quiet = quiet
        # biomarker qc
        self.biomarker_outlier_statistics = biomarker_outlier_statistics
        self.min_biomarker_threshold = min_biomarker_threshold
        self.min_number_biomarkers_for_outlier_detection = min_number_biomarkers_for_outlier_detection
        self.max_number_biomarkers_for_outlier_detection = max_number_biomarkers_for_outlier_detection
        self.report_number_of_biomarker_outlier = report_number_of_biomarker_outlier
        # Surrogate virus
        self.periode_month_surrogatevirus = periode_month_surrogatevirus
        # Water quality
        self.water_quality_number_of_last_month = water_quality_number_of_last_month
        self.min_number_of_last_measurements_for_water_qc = min_number_of_last_measurements_for_water_qc
        self.water_qc_outlier_statistics = water_qc_outlier_statistics

        self.logger = utils.SewageLogger(verbosity, quiet)
        self.__setup()

    def __setup(self):
        if not os.path.exists(self.output_folder):
            os.makedirs(self.output_folder)
        self.biomarkerQC = BiomarkerQC(self.output_folder, self.biomarker_outlier_statistics, self.min_biomarker_threshold,
                                  self.min_number_biomarkers_for_outlier_detection,
                                  self.max_number_biomarkers_for_outlier_detection,
                                  self.report_number_of_biomarker_outlier)
        self.water_quality = WaterQuality(self.output_folder, self.water_quality_number_of_last_month,
                                          self.min_number_of_last_measurements_for_water_qc, self.water_qc_outlier_statistics)

    def __check_comments(self, sample_location, measurements_df: pd.DataFrame):
        """ Flags a sewage sample if any text is stored in column 'bem_lab' or 'bem_pn' """
        flag_value = SewageFlag.COMMENT_NOT_EMPTY.value
        measurements_df[Columns.FLAG.value] += \
            np.where(((measurements_df[Columns.COMMENT_ANALYSIS.value].notnull()) & (measurements_df[Columns.COMMENT_ANALYSIS.value] != "")) |
                     ((measurements_df[Columns.COMMENT_OPERATION.value].notnull()) & (measurements_df[Columns.COMMENT_OPERATION.value] != "")),
                     flag_value,
                     0)
        num_flagged_rows = len(measurements_df[SewageFlag.is_flag_set_for_series(measurements_df[Columns.FLAG.value], SewageFlag.COMMENT_NOT_EMPTY)])
        self.logger.log.info("[Check comments] - [Sample location: '{}'] - {}/{} "
                         "samples were flagged due to non empty comments.".
                         format(sample_location, num_flagged_rows, measurements_df.shape[0]))

    def __initalize_flags(self, measurements: pd.DataFrame):
        measurements[Columns.FLAG.value] = 0
        for biomarker in Columns.get_biomarker_columns():
            measurements[Columns.get_biomarker_flag(biomarker)] = 0
        for biomarker1, biomarker2 in itertools.combinations(Columns.get_biomarker_columns(), 2):
            measurements[Columns.get_biomaker_ratio_flag(biomarker1, biomarker2)] = 0

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
        for idx, (sample_location, measurements) in enumerate(self.sewage_samples.items()):
            if idx < 3:
                continue   # skip first sewage location for testing  #Todo: remove before production
            self.logger.log.info("\n####################################################\n"
                             "\tSewage location: {} \n"
                             "####################################################".format(sample_location))
            measurements = utils.convert_sample_list2pandas(measurements)
            measurements[Columns.DATE.value] = pd.to_datetime(measurements[Columns.DATE.value],
                                                              format="%Y-%m-%d").dt.normalize()
            # Sort by collection date. Newest last.
            measurements.sort_values(by=Columns.DATE.value, ascending=True, inplace=True, ignore_index=True)
            self.logger.log.info("Total number of measurements:\t{}".format(measurements.shape[0]))

            self.__initalize_flags(measurements)

            # -----------------  BIOMARKER QC -----------------------
            # 1. check for comments. Flag samples that contain any commentary.
            self.__check_comments(sample_location, measurements)
            # 2. Mark biomarker values below threshold which are excluded from the analysis.
            self.biomarkerQC.biomarker_below_threshold(sample_location, measurements)
            # 3. Calculate pairwise biomarker values if biomarkers were not marked to be below threshold.
            self.biomarkerQC.calculate_biomarker_ratios(sample_location, measurements)
            # 4. Detect outliers
            self.biomarkerQC.detect_outliers(sample_location, measurements)
            # 5. Select measurements with enough biomarker values which are not outliers
            self.biomarkerQC.filter_required_biomarkers(sample_location, measurements)
            # 6. Create report in case the last two biomarkers were identified as outliers
            self.biomarkerQC.report_last_biomarkers_invalid(sample_location, measurements)

            # --------------------  WATER QUALITY -------------------
            self.water_quality.check_water_quality(sample_location, measurements)

            # --------------------  SUROGATVIRUS QC -------------------

            # Experimental: Final step explain flags
            measurements['flags_explained'] = SewageFlag.explain_flag_series(measurements[Columns.FLAG.value])



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
    biomarker_qc_group.add_argument('--biomarker_outlier_statistics', metavar="METHOD", default=['lof','rf','iqr'], nargs='+',
                        help=("Which outlier detection methods should be used? Multiple selections allowed. (default: 'lof','rf','iqr')\n"
                              "Possible choices are : [lof, rf, iqr, zscore, ci, all]\n"
                              "E.g. to select 'rf' and 'iqr' use: --outlier_statistics rf iqr \n"
                              "\tlof = local outlier factor\n"
                              "\trf = random forest\n"
                              "\tiqr = interquartile range\n"
                              "\tzscore = modified z-score\n"
                              "\tci = 99% confidence interval\n"
                              "\tall = use all methods\n") ,
                        choices=["lof", "rf", "iqr", "zscore", "ci", "all"],
                        required=False)
    biomarker_qc_group.add_argument('--biomarker_min_threshold', metavar="FLOAT", default=4, type=float,
                        help="Minimal biomarker threshold. (default: 4)",
                        required=False)
    biomarker_qc_group.add_argument('--min_number_biomarkers_for_outlier_detection', metavar="INT", default=9, type=int,
                        help="Minimal number of previous measurements required for outlier detection, otherwise this step is skipped. (default: 9)",
                        required=False)
    biomarker_qc_group.add_argument('--max_number_biomarkers_for_outlier_detection', metavar="INT", default=50, type=int,
                        help="Maximal number of previous measurements to use for outlier detection. (default: 50)",
                        required=False)
    biomarker_qc_group.add_argument('--report_number_of_biomarker_outliers', metavar="INT", default=2, type=int,
                        help="The number of outliers identified in the last N consecutive biomarker ratios that trigger a report. (default: 2)",
                        required=False)

    surrogatevirus_group = parser.add_argument_group("Surrogate virus quality control")
    surrogatevirus_group.add_argument('--periode_month_surrogatevirus', metavar="INT", default=4, type=int,
                        help="The periode of time (month) taken into account for surrogatevirus outliers. (default: 4)",
                        required=False)

    water_quality_group = parser.add_argument_group("Water quality")
    water_quality_group.add_argument('--water_quality_number_of_last_month', metavar="INT", default=4, type=int,
                                      help="The number of last months to be used for water quality testing. (default: 4)",
                                      required=False)
    water_quality_group.add_argument('--min_number_of_last_measurements_for_water_qc', metavar="INT", default=9, type=int,
                                     help="The minimal number of last measurements required for water quality quality control. (default: 9)",
                                     required=False)
    water_quality_group.add_argument('--water_qc_outlier_statistics', metavar="METHOD", default=['ci'], nargs='+',
                                    help=(
                                        "Which outlier detection methods should be used for water qc? Multiple selections allowed. (default: 'rf','iqr')\n"
                                        "Possible choices are : [lof, rf, iqr, zscore, ci, all]\n"
                                        "E.g. to select 'rf' and 'iqr' use: --outlier_statistics rf iqr \n"
                                        "\tlof = local outlier factor\n"
                                        "\trf = random forest\n"
                                        "\tiqr = interquartile range\n"
                                        "\tzscore = modified z-score\n"
                                        "\tci = 99% confidence interval\n"
                                        "\tall = use all methods\n"),
                                    choices=["lof", "rf", "iqr", "zscore", "ci", "all"],
                                    required=False)
    args = parser.parse_args()
    sewageQuality = SewageQuality(args.output_folder, args.verbosity, args.quiet, args.biomarker_outlier_statistics, args.biomarker_min_threshold,
                                  args.min_number_biomarkers_for_outlier_detection,
                                  args.max_number_biomarkers_for_outlier_detection,
                                  args.report_number_of_biomarker_outliers, args.periode_month_surrogatevirus,
                                  args.water_quality_number_of_last_month,
                                  args.min_number_of_last_measurements_for_water_qc, args.water_qc_outlier_statistics)
    sewageQuality.run_quality_control()

