# Created by alex at 22.06.23
import math
import numpy as np
import pandas as pd
from .utils import *
from .statistics import *
from .plotting import *


class BiomarkerQC:

    def __init__(self, output_folder, sewageStat: SewageStat, biomarker_outlier_statistics, min_biomarker_threshold,
                 min_number_biomarkers_for_outlier_detection,
                 max_number_biomarkers_for_outlier_detection, report_number_of_biomarker_outlier):
        self.output_folder = output_folder
        self.sewageStat = sewageStat
        self.biomarker_outlier_statistics = biomarker_outlier_statistics
        self.min_biomarker_threshold = min_biomarker_threshold
        self.min_number_biomarkers_for_outlier_detection = min_number_biomarkers_for_outlier_detection
        self.max_number_biomarkers_for_outlier_detection = max_number_biomarkers_for_outlier_detection
        self.report_number_of_biomarker_outlier = report_number_of_biomarker_outlier
        self.logger = SewageLogger(self.output_folder)


    def check_comments(self, sample_location: str, measurements: pd.DataFrame, index):
        """ Flags a sewage sample if any text is stored in column 'bem_lab' or 'bem_pn' """
        current_measurement = measurements.iloc[index]
        comment_analysis = current_measurement[Columns.COMMENT_ANALYSIS.value]
        comment_operation = current_measurement[Columns.COMMENT_OPERATION.value]
        if comment_analysis != '' or comment_operation != '':
            SewageFlag.add_flag_to_index_column(measurements, index, CalculatedColumns.FLAG.value, SewageFlag.COMMENT_NOT_EMPTY)
            self.sewageStat.add_comment_not_empty()

    def check_mean_sewage_flow_present(self, sample_location, measurements: pd.DataFrame, index):
        current_measurement = measurements.iloc[index]
        mean_sewage_flow = current_measurement[Columns.MEAN_SEWAGE_FLOW.value]
        if math.isnan(mean_sewage_flow) or mean_sewage_flow == 0:
            SewageFlag.add_flag_to_index_column(measurements, index, CalculatedColumns.FLAG.value, SewageFlag.MISSING_MEAN_SEWAGE_FLOW)
            self.sewageStat.add_mean_sewage_flow_empty()


    def biomarker_below_threshold_or_empty(self, sample_location, measurements: pd.DataFrame, index):
        """
        Mark biomarker values that are below given threshold or empty.
        """
        current_measurement = measurements.iloc[index]
        for biomarker in Columns.get_biomarker_columns():
            biomarker_value = current_measurement[biomarker]
            if math.isnan(biomarker_value) or biomarker_value < self.min_biomarker_threshold:
                SewageFlag.add_flag_to_index_column(measurements, index, CalculatedColumns.get_biomarker_flag(biomarker), SewageFlag.BIOMARKER_BELOW_THRESHOLD_OR_EMPTY)
                self.sewageStat.add_biomarker_below_threshold_or_empty(biomarker)

    def standardize_biomarker_values(self, sample_location,  measurements: pd.DataFrame):
        for biomarker in Columns.get_biomarker_columns():
            index = measurements.columns.get_loc(biomarker)
            standardized_biomarker = ( measurements[biomarker] - measurements[biomarker].mean()) / measurements[biomarker].std()
            measurements.insert(index + 1, biomarker+"_zscore", standardized_biomarker)

    def calculate_biomarker_ratios(self, sample_location, measurements: pd.DataFrame, index):
        """
        Calculate pairwise biomarker ratios only if both are above the minimal threshold level.
        In case one of the biomarker values is below the threshold, the ratio is not calculated.
        """
        current_measurement = measurements.iloc[index]
        for biomarker1, biomarker2 in itertools.combinations(Columns.get_biomarker_columns(), 2):

            is_biomarker1_flagged = SewageFlag.is_flag(current_measurement[CalculatedColumns.get_biomarker_flag(biomarker1)],
                                                       SewageFlag.BIOMARKER_BELOW_THRESHOLD_OR_EMPTY)
            is_biomarker2_flagged = SewageFlag.is_flag(current_measurement[CalculatedColumns.get_biomarker_flag(biomarker2)],
                                                       SewageFlag.BIOMARKER_BELOW_THRESHOLD_OR_EMPTY)
            biomarker1_value, biomarker2_value = current_measurement[biomarker1 + "_zscore"], current_measurement[biomarker2 + "_zscore"]
            if is_biomarker1_flagged or is_biomarker2_flagged or biomarker1_value == 0 or biomarker2_value == 0:
                measurements.at[index, biomarker1 + "/" + biomarker2] = np.NAN
            else:
                last_biomarker_ratios = self.__get_previous_biomarkers_ratios(measurements, index, biomarker1, biomarker2)
                last_biomarker_ratio_median = np.median(last_biomarker_ratios[biomarker1+"_zscore"] / last_biomarker_ratios[biomarker2+"_zscore"]) if len(last_biomarker_ratios) > 0 else np.NAN
                if not np.isnan(last_biomarker_ratio_median):
                    measurements.at[index, biomarker1 + "/" + biomarker2] = (biomarker1_value / biomarker2_value) / last_biomarker_ratio_median
                else:
                    measurements.at[index, biomarker1 + "/" + biomarker2] = 1


    def __get_previous_biomarkers_ratios(self, measurements_df: pd.DataFrame, index, biomarker1, biomarker2):
        """
          Obtain last biomarker ratios for the previous measurements starting from current measurement.
          Remove previously detected outliers and empty ratios.
        """
        max_idx = (index - 1 - self.max_number_biomarkers_for_outlier_detection) \
            if (index - 1 - self.max_number_biomarkers_for_outlier_detection) >= 0 else 0
        min_idx = index if index >= 0 else 0
        last_N_measurements = measurements_df.iloc[max_idx: min_idx]
        # remove previously detected outliers
        biomarker_ratio_flags = last_N_measurements[CalculatedColumns.get_biomaker_ratio_flag(biomarker1, biomarker2)]
        last_N_measurements = last_N_measurements[SewageFlag.is_not_flag_set_for_series(biomarker_ratio_flags,
                                                                                        SewageFlag.BIOMARKER_RATIO_OUTLIER)]
        biomarker_ratios = last_N_measurements[[Columns.DATE.value, biomarker1 + "/" + biomarker2, biomarker1+"_zscore", biomarker2+"_zscore"]]
        # remove empty ratios
        biomarker_ratios = biomarker_ratios.dropna()
        return biomarker_ratios

    def detect_outliers(self, sample_location, measurements: pd.DataFrame, index):
        """
           For each measurement and biomarker pair outliers will be marked.
        """
        current_measurement = measurements.iloc[index]
        for biomarker1, biomarker2 in itertools.combinations(Columns.get_biomarker_columns(), 2):
            biomarker_ratio = biomarker1 + "/" + biomarker2
            if current_measurement[biomarker_ratio] and not math.isnan(current_measurement[biomarker_ratio]):  # only if biomarker ratio is not None for current measurement
                last_biomarker_ratios = self.__get_previous_biomarkers_ratios(measurements, index, biomarker1, biomarker2)
                if len(last_biomarker_ratios) < self.min_number_biomarkers_for_outlier_detection:
                    SewageFlag.add_flag_to_index_column(measurements, index, CalculatedColumns.get_biomaker_ratio_flag(biomarker1, biomarker2),
                                                        SewageFlag.NOT_ENOUGH_PREVIOUS_BIOMARKER_VALUES)
                    continue
                is_outlier = detect_outliers(self.biomarker_outlier_statistics, last_biomarker_ratios[biomarker_ratio], current_measurement[biomarker_ratio], isFactor=True)
                if is_outlier:
                    SewageFlag.add_flag_to_index_column(measurements, index, CalculatedColumns.get_biomaker_ratio_flag(biomarker1, biomarker2), SewageFlag.BIOMARKER_RATIO_OUTLIER)
                    self.sewageStat.add_biomarker_ratio_outlier(biomarker1, biomarker2, 'outlier')
                else:
                    self.sewageStat.add_biomarker_ratio_outlier(biomarker1, biomarker2, 'passed')
            else:
                self.sewageStat.add_biomarker_ratio_outlier(biomarker1, biomarker2, 'skipped')


    def __get_usable_biomarkers(self, row: pd.Series) -> ([], []):
        """
        filters biomarkers which are usable, thus either empty or below threshold
        """
        row = row[CalculatedColumns.get_biomarker_flag_columns()]
        usable_biomarkers = SewageFlag.is_flag_set_for_series(row, SewageFlag.BIOMARKER_BELOW_THRESHOLD_OR_EMPTY)
        usable_biomarkers = usable_biomarkers[usable_biomarkers == False]
        usable_biomarkers = list(map(lambda x: x.replace(CalculatedColumns.BIOMARKER_FLAG.value + "_", ""), list(usable_biomarkers.index)))
        biomarker_ratio_flags = []
        for biomarker1, biomarker2 in itertools.combinations(Columns.get_biomarker_columns(), 2):
            if biomarker1 in usable_biomarkers and biomarker2 in usable_biomarkers:
                biomarker_ratio_flags.append(CalculatedColumns.get_biomaker_ratio_flag(biomarker1, biomarker2))
        return usable_biomarkers, biomarker_ratio_flags

    def __mark_biomakers_with_all_outliers(self, current_measurement):
        usable_biomarkers, biomarker_ratio_flags = self.__get_usable_biomarkers(current_measurement)
        biomarker_outlier_set = set()
        if len(usable_biomarkers) > 2:
            for biomarker in usable_biomarkers:
                biomarker_ratio_flag_df = current_measurement[[b for b in biomarker_ratio_flags if biomarker in b]]
                outliers_for_biomarker = SewageFlag.is_flag_set_for_series(biomarker_ratio_flag_df, SewageFlag.BIOMARKER_RATIO_OUTLIER)
                if outliers_for_biomarker.all():
                    biomarker_outlier_set.add(biomarker)
        return biomarker_outlier_set

    def __flag_outliers_with_all_ratio_outliers(self, biomarker_outlier_set: set, index, measurements_df: pd.DataFrame) -> None:
        for biomarker in biomarker_outlier_set:
            SewageFlag.add_flag_to_index_column(measurements_df, index, CalculatedColumns.get_biomarker_flag(biomarker),
                                                SewageFlag.BIOMARKER_VALIDATED_OUTLIER)

    def assign_biomarker_outliers_based_on_ratio_flags(self, sample_location, measurements: pd.DataFrame, index):
        current_measurement = measurements.iloc[index]
        usable_biomarkers, biomarker_ratio_flags = self.__get_usable_biomarkers(current_measurement)
        all_ratios_outlier_biomarker_set = self.__mark_biomakers_with_all_outliers(current_measurement)
        # if all ratios for a given biomarker are flagged as outliers then the biomarkers is flagged as validated outlier
        self.__flag_outliers_with_all_ratio_outliers(all_ratios_outlier_biomarker_set, index, measurements)
        for biomarker1, biomarker2 in itertools.combinations(usable_biomarkers, 2):
            # is biomarker1 or biomarker2 is not already flagged as outlier before
            if not biomarker1 in all_ratios_outlier_biomarker_set and not biomarker2 in all_ratios_outlier_biomarker_set:
                is_outlier = SewageFlag.is_flag(current_measurement[CalculatedColumns.get_biomaker_ratio_flag(biomarker1, biomarker2)], SewageFlag.BIOMARKER_RATIO_OUTLIER)
                if is_outlier:
                    SewageFlag.add_flag_to_index_column(measurements, index, CalculatedColumns.get_biomarker_flag(biomarker1), SewageFlag.BIOMARKER_PROBABLE_OUTLIER)
                    SewageFlag.add_flag_to_index_column(measurements, index, CalculatedColumns.get_biomarker_flag(biomarker2), SewageFlag.BIOMARKER_PROBABLE_OUTLIER)

    def analyze_usable_biomarkers(self, sample_location, measurements: pd.DataFrame, index):
        """
        Sums the number of biomarkers which are not marked as outliers using the outlier columns.
        Flags measurments which have less than specified number of valid biomarkers.
        """
        current_measurement = measurements.iloc[index]
        num_usable_biomarkers = 0
        is_probable_outlier_flag = False
        for biomarker in Columns.get_biomarker_columns():
            is_below_threshold_or_empty = SewageFlag.is_flag(current_measurement[CalculatedColumns.get_biomarker_flag(biomarker)], SewageFlag.BIOMARKER_BELOW_THRESHOLD_OR_EMPTY)
            is_validated_outlier = SewageFlag.is_flag(current_measurement[CalculatedColumns.get_biomarker_flag(biomarker)], SewageFlag.BIOMARKER_VALIDATED_OUTLIER)
            is_probable_outlier = SewageFlag.is_flag(current_measurement[CalculatedColumns.get_biomarker_flag(biomarker)], SewageFlag.BIOMARKER_PROBABLE_OUTLIER)
            if not is_below_threshold_or_empty and not is_validated_outlier: #and not is_probable_outlier:
                num_usable_biomarkers += 1
            if is_probable_outlier:
                is_probable_outlier_flag = True
        measurements.at[index, CalculatedColumns.NUMBER_OF_USABLE_BIOMARKERS.value] = num_usable_biomarkers
        if is_probable_outlier_flag:  # in case any biomarker was marked as probable outlier set flag for normalization
            SewageFlag.add_flag_to_index_column(measurements, index, CalculatedColumns.FLAG.value, SewageFlag.BIOMARKER_PROBABLE_OUTLIER)


    def report_last_biomarkers_invalid(self, sample_location, measurements_df: pd.DataFrame):
        last_two_measurements = measurements_df.tail(self.report_number_of_biomarker_outlier)
        if last_two_measurements.shape[0] == self.report_number_of_biomarker_outlier:
            for biomarker1, biomarker2 in itertools.combinations(Columns.get_biomarker_columns(), 2):
                valid_ratios = SewageFlag.is_not_flag_set_for_series(last_two_measurements[CalculatedColumns.get_biomaker_ratio_flag(biomarker1, biomarker2)],
                                                                     SewageFlag.BIOMARKER_RATIO_OUTLIER)
                num_valid_ratios_last_measurements = valid_ratios.sum()
                if num_valid_ratios_last_measurements == 0:  # last samples where marked as outlier
                    biomarker_last_outliers = last_two_measurements[
                        [Columns.DATE.value, biomarker1, biomarker2, biomarker1 + "/" + biomarker2]]
                    self.logger.log.warn("[Report last two biomarker outliers] - [Sample location: '{}'] - "
                                         "The last two measurements for biomarker ratio '{}/{}' are marked as outliers.\n"
                                         "{}".format(sample_location, biomarker1, biomarker2, biomarker_last_outliers))
                    # Todo: handle reporting
        else:
            self.logger.log.debug("[Report last two biomarker outliers] - [Sample location: '{}'] - "
                                  "Less than {} last values available. Skipping... ".format(sample_location, self.report_number_of_biomarker_outlier))
