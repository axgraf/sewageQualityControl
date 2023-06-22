# Created by alex at 22.06.23
import itertools
import numpy as np
from .flags import Flag
from .utils import *

class BiomarkerQC:

    def __init__(self, min_biomarker_threshold, minimal_number_measurements_for_outlier_detection,
                 max_number_measurements_for_outlier_detection, min_number_of_valid_biomarkers):
        self.min_biomarker_threshold = min_biomarker_threshold
        self.minimal_number_measurements_for_outlier_detection = minimal_number_measurements_for_outlier_detection
        self.max_number_measurements_for_outlier_detection = max_number_measurements_for_outlier_detection
        self.min_number_of_valid_biomarkers = min_number_of_valid_biomarkers
        self.biomarker_columns = ['biomarker_N1', 'biomarker_N2', 'biomarker_N3', 'biomarker_E', 'biomarker_ORF',
                                  'biomarker_RDRP']
        self.logger = SewageLogger()


    def biomarker_below_threshold(self, sample_location, measurements_df: pd.DataFrame):
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
        self.logger.log.debug("[Biomarker below threshold] - [Sample location: '{}'] - \n"
                         "\tNumber of excluded biomarkers which are below threshold or empty:\n\t\t{}\n"
                         "\tout of {} samples.".
                         format(sample_location, biomarkers_below_dict, measurements_df.shape[0]))

    def calculate_biomarker_ratios(self, sample_location, measurements_df: pd.DataFrame):
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

    def __get_previous_biomarkers(self, measurements_df: pd.DataFrame, index, biomarker1, biomarker2):
        """
          Obtain last biomarker ratios for the previous measurements starting from current measurement.
          Remove previously detected outliers and empty ratios.
        """
        max_idx = (index - 1 - self.max_number_measurements_for_outlier_detection) \
            if (index - 1 - self.max_number_measurements_for_outlier_detection) >= 0 else 0
        min_idx = index if index >= 0 else 0
        last_N_measurements = measurements_df.iloc[max_idx: min_idx]
        if biomarker1 + "/" + biomarker2 + "_outlier" in last_N_measurements:  # remove previously detected outliers
            last_N_measurements = last_N_measurements[
                last_N_measurements[biomarker1 + "/" + biomarker2 + "_outlier"] != True]
        else:
            self.logger.log.debug("No outliers detected yet")
        biomarker_ratios = last_N_measurements[biomarker1 + "/" + biomarker2]
        # remove empty ratios
        biomarker_ratios = biomarker_ratios.dropna()
        return biomarker_ratios

    def __is_biomarker_outlier(self, sample_location, collectionDate, current_biomarker_ratio, last_biomarker_ratios, biomarker1, biomarker2):
        """
            Detect outliers for a given biomarker set using previous measurements with ...
        """
        # Todo: Which outlier detection should we use?
        is_confidence_interval_99_outlier, confidence_interval = is_confidence_interval_outlier(
            current_biomarker_ratio, last_biomarker_ratios, 0.99)
        is_iqr_outlier, iqr_range = inter_quantil_range(current_biomarker_ratio, last_biomarker_ratios)
        is_zscore_outlier, z_score_threshold = is_outlier_modified_z_score(current_biomarker_ratio,
                                                                                 last_biomarker_ratios)
        self.logger.log.debug("[Biomarker outlier detection] - [Sample location: '{}'] - [Collection date: '{}']\n"
                          "\tBiomarker ratio: {}\t\n"
                          "\tCurrent ratio:\t{}\n"
                          "\tLast ratios:\t{}\n"
                          "\tIs Z-score-outlier:\t\t{}\t(Threshold: {} < 3.5)\n"
                          "\tIs IQR-outlier:\t{}\t\tRange: {}\n"
                          "\tIs 99% confidence interval outlier:\t{}\tRange: {}".format(sample_location,
                                collectionDate, biomarker1 + "/" + biomarker2, current_biomarker_ratio,
                                ['%.4f' % elem for elem in last_biomarker_ratios.tolist()], is_zscore_outlier, z_score_threshold,
                                is_iqr_outlier, iqr_range, is_confidence_interval_99_outlier, confidence_interval))
        if is_zscore_outlier and is_iqr_outlier:
            return True
        return False

    def detect_outliers(self, sample_location, measurements_df: pd.DataFrame):
        """
           For each measurement and biomarker pair outliers will be marked.
        """
        skipped = 0
        stat_dict = dict()
        for index, current_measurement in measurements_df.iterrows():
            for biomarker1, biomarker2 in itertools.combinations(self.biomarker_columns, 2):
                add_default_biomarker_statistic(stat_dict, biomarker1, biomarker2)
                current_biomarker_ratio = current_measurement[biomarker1 + "/" + biomarker2]
                if current_biomarker_ratio:  # only if biomarker ratio is not None for current measurement
                    stat_dict[biomarker1 + "/" + biomarker2]["total"] += 1
                    last_biomarker_ratios = self.__get_previous_biomarkers(measurements_df, index, biomarker1, biomarker2)
                    measurements_df.at[index, Flag.NOT_ENOUGH_BIOMARKERS_FOR_OUTLIER_DETECTION.name] = False
                    if len(last_biomarker_ratios) < self.minimal_number_measurements_for_outlier_detection:
                        measurements_df.at[index, Flag.NOT_ENOUGH_BIOMARKERS_FOR_OUTLIER_DETECTION.name] = True
                        self.logger.log.debug("[Biomarker outlier detection] - [Sample location: '{}'] - [Collection date: '{}'] - "
                                            "Skipping biomarker outlier detection. "
                                            "\n\tFound '{}' previous measurements. "
                                            "More than '{}' previous measurements are required.".format(sample_location, current_measurement['collectionDate'],
                                                                                                        len(last_biomarker_ratios),
                                                                                                        self.minimal_number_measurements_for_outlier_detection))
                        stat_dict[biomarker1 + "/" + biomarker2]["skipped"] += 1
                        continue
                    is_outlier = self.__is_biomarker_outlier(sample_location, current_measurement['collectionDate'],
                                                             current_biomarker_ratio, last_biomarker_ratios,
                                                             biomarker1, biomarker2)
                    if is_outlier:
                        stat_dict[biomarker1 + "/" + biomarker2]["failed"] += 1
                        measurements_df.at[index, biomarker1 + "/" + biomarker2 + "_outlier"] = True
                    else:
                        stat_dict[biomarker1 + "/" + biomarker2]["passed"] += 1
                        measurements_df.at[index, biomarker1 + "/" + biomarker2 + "_outlier"] = False
                else:
                    self.logger.log.debug("[Biomarker outlier detection] - [Sample location: '{}'] - [Collection date: '{}'] - "
                                      "Empty biomarker ratio '{}/{}' with value: '{}'. Biomarker ratio will be skipped."
                                      .format(sample_location, current_measurement['collectionDate'],
                                                                                biomarker1, biomarker2,
                                                                                current_biomarker_ratio))
        self.logger.log.info("[Biomarker outlier detection] - [Sample location: '{}'] - \n"
                         "Biomarkers ratios:\n{}".format(sample_location, pretty_print_biomarker_statistic(stat_dict)))


    def select_minimal_required_biomarkers(self, sample_location, measurements_df: pd.DataFrame):
        """
        Sums the number of biomarkers which are not marked as outliers using the outlier columns.
        Flags measurments which have less than specified number of valid biomarkers.
        """
        outlier_columns = []
        for biomarker1, biomarker2 in itertools.combinations(self.biomarker_columns, 2):
            outlier_columns.append(biomarker1 + "/" + biomarker2 + "_outlier")
        outlier_frame = measurements_df.loc[:, [i in outlier_columns for i in measurements_df.columns]]
        measurements_df['biomarkerSums'] = (outlier_frame == False).sum(axis=1)
        measurements_df[Flag.MIN_BIOMARKER_NUMBER_NOT_REACHED.name] = \
            np.where(measurements_df['biomarkerSums'] >= self.min_number_of_valid_biomarkers, True, False)
        count_biomarker_ready_for_analysis = measurements_df[Flag.MIN_BIOMARKER_NUMBER_NOT_REACHED.name].value_counts()
        self.logger.log.info("[Select Biomarkers] - [Sample location: '{}'] - \n"
                         "Biomarkers with >= '{}' valid values:\n{}".format(sample_location,
                                                                            self.min_number_of_valid_biomarkers,
                                                                            count_biomarker_ready_for_analysis.to_string()))


