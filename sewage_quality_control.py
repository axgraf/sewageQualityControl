import itertools
import argparse
import numpy as np
import datetime
from lib.arcgis import *
from lib.flags import Flag
import lib.utils as utils
import logging
import scipy.stats as st


# logging.basicConfig(format='%(asctime)s - %(levelname)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S',
#                    stream=sys.stdout, level=logging.INFO)


class SewageQuality:

    def __init__(self, verbosity, quiet, min_biomarker_threshold, minimal_number_measurements_for_outlier_detection,
                 max_number_measurements_for_outlier_detection):
        self.sewage_samples = None
        self.verbosity = verbosity
        self.quiet = quiet
        self.min_biomarker_threshold = min_biomarker_threshold
        self.minimal_number_measurements_for_outlier_detection = minimal_number_measurements_for_outlier_detection
        self.max_number_measurements_for_outlier_detection = max_number_measurements_for_outlier_detection
        self.biomarker_columns = ['biomarker_N1', 'biomarker_N2', 'biomarker_N3', 'biomarker_E', 'biomarker_ORF',
                                  'biomarker_RDRP']
        self.__create_logger()

    def __create_logger(self):
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.WARNING)
        if self.verbosity and self.verbosity == 1:
            self.logger.setLevel(logging.INFO)
        elif self.verbosity and self.verbosity >= 2:
            self.logger.setLevel(logging.DEBUG)
        if self.quiet:
            self.logger.setLevel(logging.ERROR)
        fmt = "%(asctime)s - %(levelname)s - %(message)s"
        stdout_handler = logging.StreamHandler()
        stdout_handler.setLevel(logging.DEBUG)
        stdout_handler.setFormatter(utils.CustomFormatter(fmt))
        today = datetime.date.today()
        file_handler = logging.FileHandler('sewage_quality_{}.log'.format(today.strftime('%Y_%m_%d')))
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(logging.Formatter(fmt))
        self.logger.addHandler(stdout_handler)
        self.logger.addHandler(file_handler)

    def __check_comments(self, sample_location, measurements_df):
        """ Flags a sewage sample if any not was given in column 'bem_lab' or 'bem_pn' """
        measurements_df[Flag.COMMENT_NOT_EMPTY.name] = \
            np.where(((measurements_df['bem_lab'].notnull()) & (measurements_df['bem_lab'] != "")) |
                     ((measurements_df['bem_pn'].notnull()) & (measurements_df['bem_pn'] != "")), True, False)
        num_flagged_rows = len(measurements_df[measurements_df[Flag.COMMENT_NOT_EMPTY.name] == True])
        self.logger.info("[Check comments] - [Sample location: '{}'] - {}/{} "
                         "samples were flagged due to non empty comments.".
                         format(sample_location, num_flagged_rows, measurements_df.shape[0]))

    def __biomarker_below_threshold(self, sample_location, measurements_df):
        """
        Mark biomarker values that are below given threshold or empty.
        """
        biomarkers_below_dict = dict()
        for biomarker in self.biomarker_columns:
            measurements_df[biomarker + "_below_threshold"] = np.where(
                ((measurements_df[biomarker].isnull()) | (measurements_df[biomarker] < self.min_biomarker_threshold)),
                True, False)
            biomarkers_below_dict[biomarker] = len(
                measurements_df[measurements_df[biomarker + "_below_threshold"] == True])
        self.logger.debug("[Biomarker below threshold] - [Sample location: '{}'] - \n"
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
                np.where(
                    measurements_df[biomarker1 + "_below_threshold"] | measurements_df[biomarker2 + "_below_threshold"],
                    None,  # if one biomarker is below threshold
                    measurements_df[biomarker1] / measurements_df[biomarker2])

    def __get_last_biomarkers_for_last_N_measurements(self, measurements_df, index, biomarker1, biomarker2):
        max_idx = (index - 1 - self.max_number_measurements_for_outlier_detection) \
            if (index - 1 - self.max_number_measurements_for_outlier_detection) >= 0 else 0
        min_idx = index if index >= 0 else 0
        last_N_measurements = measurements_df.iloc[max_idx: min_idx]
        if biomarker1 + "/" + biomarker2 + "_outlier" in last_N_measurements:  # remove previously detected outliers
            last_N_measurements = last_N_measurements[
                last_N_measurements[biomarker1 + "/" + biomarker2 + "_outlier"] != True]
        else:
            self.logger.debug("no confidence interval calculated yet")
        biomarker_ratios = last_N_measurements[biomarker1 + "/" + biomarker2]
        # remove empty ratios
        biomarker_ratios = biomarker_ratios.dropna()
        return biomarker_ratios

    def __detect_outliers(self, sample_location, measurements_df):
        for index, current_measurement in measurements_df.iterrows():
            for biomarker1, biomarker2 in itertools.combinations(self.biomarker_columns, 2):
                current_biomarker_ratio = current_measurement[biomarker1 + "/" + biomarker2]
                if current_biomarker_ratio:  # only if biomarker ratio is not None for current measurement
                    biomarker_ratios = self.__get_last_biomarkers_for_last_N_measurements(
                        measurements_df, index, biomarker1, biomarker2)
                    measurements_df.at[index, Flag.NOT_ENOUGH_BIOMARKERS_FOR_OUTLIER_DETECTION.name] = False
                    if len(biomarker_ratios) < self.minimal_number_measurements_for_outlier_detection:
                        measurements_df.at[index, Flag.NOT_ENOUGH_BIOMARKERS_FOR_OUTLIER_DETECTION.name] = True
                        self.logger.debug("[Biomarker outlier detection] - [Sample location: '{}'] - [Collection date: '{}'] - "
                                            "Skipping biomarker outlier detection. "
                                            "\n\tFound '{}' previous measurements. "
                                            "More than '{}' previous measurements are required.".format(sample_location, current_measurement['collectionDate'],
                                                                                                        len(biomarker_ratios),
                                                                                                        self.minimal_number_measurements_for_outlier_detection))
                        continue
                    # Todo: Which outlier detection should we use?
                    is_confidence_interval_99_outlier, confidence_interval = utils.is_confidence_interval_outlier(
                        current_biomarker_ratio, biomarker_ratios, 0.99)
                    is_iqr_outlier, iqr_range = utils.inter_quantil_range(current_biomarker_ratio, biomarker_ratios)
                    is_zscore_outlier, z_score_threshold = utils.is_outlier_modified_z_score(current_biomarker_ratio,
                                                                                   biomarker_ratios)
                    self.logger.debug("[Biomarker outlier detection] - [Sample location: '{}'] - [Collection date: '{}']\n"
                                      "\tBiomarker ratio: {}\t\n"
                                      "\tCurrent ratio:\t{}\n"
                                      "\tLast ratios:\t{}\n"
                                      "\tIs Z-score-outlier:\t\t{}\t(Threshold: {} < 3.5)\n"
                                      "\tIs IQR-outlier:\t{}\t\tRange: {}\n"
                                      "\tIs 99% confidence interval outlier:\t{}\tRange: {}".format(sample_location, current_measurement['collectionDate'],
                        biomarker1+"/"+biomarker2, current_biomarker_ratio, ['%.4f' % elem for elem in biomarker_ratios.tolist()],
                        is_zscore_outlier, z_score_threshold, is_iqr_outlier, iqr_range, is_confidence_interval_99_outlier, confidence_interval))

                    if is_zscore_outlier and is_iqr_outlier:
                        measurements_df.at[index, biomarker1 + "/" + biomarker2 + "_outlier"] = True
                    else:
                        measurements_df.at[index, biomarker1 + "/" + biomarker2 + "_outlier"] = False
                else:
                    self.logger.debug("[Biomarker outlier detection] - [Sample location: '{}'] - [Collection date: '{}'] - "
                                      "Biomarker ratio '{}/{}' with value: '{}' not used as it "
                                      "not calculated in current sample".format(sample_location, current_measurement['collectionDate'],
                                                                                biomarker1, biomarker2,
                                                                                current_biomarker_ratio))

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
            self.logger.info("########################  Sewage location: {}  ########################\n".format(sample_location))
            measurements = utils.convert_sample_list2pandas(measurements)
            measurements.sort_values(by='collectionDate', ascending=True, inplace=True,
                                     ignore_index=True)  # Sort by collection date. Newest first.
            # -----------------  BIOMARKER QC -----------------------
            # 1. check for comments. Flag samples that contain any commentary
            self.__check_comments(sample_location, measurements)
            # 2. Mark biomarker values below threshold which are excluded from the analysis.
            self.__biomarker_below_threshold(sample_location, measurements)
            # 3. Calculate pairwise biomarker values if biomarkers were not marked to be below threshold.
            self.__calculate_biomarker_ratios(sample_location, measurements)
            # 4. Detect outliers
            self.__detect_outliers(sample_location, measurements)
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
    parser.add_argument('-m', '--minimal_number_measurements_for_outlier_detection', metavar="INT", default=9, type=int,
                        help="Minimal number of previous measurements required for outlier detection, otherwise this step is skipped",
                        required=False)
    parser.add_argument('-x', '--max_number_measurements_for_outlier_detection', metavar="INT", default=50, type=int,
                        help="Maximal number of previous measurements to use for outlier detection",
                        required=False)
    parser.add_argument('-v', '--verbosity', action="count", help="Increase output verbosity. Can be used multiple times.")
    parser.add_argument('-q', '--quiet', action='store_true', help="Print litte output.")
    args = parser.parse_args()
    sewageQuality = SewageQuality(args.verbosity, args.quiet, args.biomarker_min_threshold,
                                  args.minimal_number_measurements_for_outlier_detection,
                                  args.max_number_measurements_for_outlier_detection)
    sewageQuality.run_quality_control()
