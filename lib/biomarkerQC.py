# Created by alex at 22.06.23
import math
import numpy as np
import pandas as pd
from .utils import *
from .plotting import *


class BiomarkerQC:

    def __init__(self, output_folder, interactive, biomarker_outlier_statistics, min_biomarker_threshold,
                 min_number_biomarkers_for_outlier_detection,
                 max_number_biomarkers_for_outlier_detection, report_number_of_biomarker_outlier, plotting=True):
        self.output_folder = output_folder
        self.interactive = interactive
        self.biomarker_outlier_statistics = biomarker_outlier_statistics
        self.min_biomarker_threshold = min_biomarker_threshold
        self.min_number_biomarkers_for_outlier_detection = min_number_biomarkers_for_outlier_detection
        self.max_number_biomarkers_for_outlier_detection = max_number_biomarkers_for_outlier_detection
        self.report_number_of_biomarker_outlier = report_number_of_biomarker_outlier
        self.plotting = plotting
        self.logger = SewageLogger(self.output_folder)

    def check_comments(self, sample_location, measurements_df: pd.DataFrame):
        """ Flags a sewage sample if any text is stored in column 'bem_lab' or 'bem_pn' """
        flag_value = SewageFlag.COMMENT_NOT_EMPTY.value
        measurements_df[Columns.COMMENT_ANALYSIS.value] = measurements_df[Columns.COMMENT_ANALYSIS.value].astype(str).str.strip()
        measurements_df[Columns.COMMENT_OPERATION.value] = measurements_df[Columns.COMMENT_OPERATION.value].astype(str).str.strip()
        flag_series = \
            np.where(((measurements_df[Columns.COMMENT_ANALYSIS.value].notnull()) & (
                    measurements_df[Columns.COMMENT_ANALYSIS.value] != "")) |
                     ((measurements_df[Columns.COMMENT_OPERATION.value].notnull()) & (
                             measurements_df[Columns.COMMENT_OPERATION.value] != "")),
                     flag_value,
                     0)
        SewageFlag.add_series_flag_to_column(measurements_df, CalculatedColumns.FLAG.value, flag_series)
        num_flagged_rows = len(measurements_df[SewageFlag.is_flag_set_for_series(measurements_df[CalculatedColumns.FLAG.value],
                                                                                 SewageFlag.COMMENT_NOT_EMPTY)])
        self.logger.log.info("[Check comments] - [Sample location: '{}'] - {}/{} "
                             "samples were flagged due to non empty comments.".
                             format(sample_location, num_flagged_rows, measurements_df.shape[0]))

    def check_mean_sewage_flow_present(self, sample_location, measurements_df: pd.DataFrame):
        flag_series = np.where((measurements_df[Columns.MEAN_SEWAGE_FLOW.value].isna()),
                               SewageFlag.MISSING_MEAN_SEWAGE_FLOW.value,
                               0)
        SewageFlag.add_series_flag_to_column(measurements_df, CalculatedColumns.FLAG.value, flag_series)
        num_flagged_rows = len(measurements_df[SewageFlag.is_flag_set_for_series(measurements_df[CalculatedColumns.FLAG.value],
                                                                                 SewageFlag.MISSING_MEAN_SEWAGE_FLOW)])
        self.logger.log.info("[Check mean sewage flow] - [Sample location: '{}'] - {}/{} "
                             "samples were flagged due to missing mean sewage flow.".
                             format(sample_location, num_flagged_rows, measurements_df.shape[0]))

    def biomarker_below_threshold_or_empty(self, sample_location, measurements_df: pd.DataFrame):
        """
        Mark biomarker values that are below given threshold or empty.
        """
        biomarkers_below_dict = dict()
        flag_value = SewageFlag.BIOMARKER_BELOW_THRESHOLD_OR_EMPTY.value
        for biomarker in Columns.get_biomarker_columns():
            measurements_df[CalculatedColumns.get_biomarker_flag(biomarker)] += np.where(
                ((measurements_df[biomarker].isnull()) | (measurements_df[biomarker] < self.min_biomarker_threshold)),
                flag_value, 0)
            biomarker_flag_series = measurements_df[CalculatedColumns.get_biomarker_flag(biomarker)]
            biomarkers_below_dict[biomarker] = \
                len(measurements_df[SewageFlag.is_flag_set_for_series(biomarker_flag_series, SewageFlag.BIOMARKER_BELOW_THRESHOLD_OR_EMPTY)])
        self.logger.log.info("[Biomarker below threshold] - [Sample location: '{}'] - \n"
                             "\tNumber of excluded biomarkers which are below threshold ( < {} ) or empty:\n\t{}\n"
                             "\tout of {} samples.".
                             format(sample_location, self.min_biomarker_threshold, biomarkers_below_dict,
                                    measurements_df.shape[0]))

    def calculate_biomarker_ratios(self, sample_location, measurements_df: pd.DataFrame):
        """
        Calculate pairwise biomarker ratios only if both are above the minimal threshold level.
        In case one of the biomarker values is below the threshold, the ratio is not calculated.
        """
        biomarkers_ratio_dict = dict()
        for biomarker1, biomarker2 in itertools.combinations(Columns.get_biomarker_columns(), 2):
            biomarker1_series = measurements_df[CalculatedColumns.get_biomarker_flag(biomarker1)]
            biomarker2_series = measurements_df[CalculatedColumns.get_biomarker_flag(biomarker2)]
            with np.errstate(divide='ignore'):  # ignore zero divisions as they were converted to Nan
                measurements_df[biomarker1 + "/" + biomarker2] = \
                    np.where((SewageFlag.is_flag_set_for_series(biomarker1_series, SewageFlag.BIOMARKER_BELOW_THRESHOLD_OR_EMPTY)) |
                             (SewageFlag.is_flag_set_for_series(biomarker2_series, SewageFlag.BIOMARKER_BELOW_THRESHOLD_OR_EMPTY)),
                             np.NAN,  # if one biomarker is below threshold no threshold will be calculated
                             np.log2(list(measurements_df[biomarker1] / measurements_df[biomarker2])))
            num_valid_ratios = len(measurements_df[measurements_df[biomarker1 + "/" + biomarker2].notnull()])
            if num_valid_ratios > 0:
                biomarkers_ratio_dict[biomarker1 + "/" + biomarker2] = num_valid_ratios
        self.logger.log.info("[Calculate biomarker ratios] - [Sample location: '{}'] -\n"
                             "\tNumber of calculated biomarker ratios:\n"
                             "\t{}".format(sample_location, biomarkers_ratio_dict))

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
        biomarker_ratios = last_N_measurements[[Columns.DATE.value, biomarker1 + "/" + biomarker2, biomarker1, biomarker2]]
        # remove empty ratios
        biomarker_ratios = biomarker_ratios.dropna()
        return biomarker_ratios

    def detect_outliers(self, sample_location, measurements_df: pd.DataFrame):
        """
           For each measurement and biomarker pair outliers will be marked.
        """
        self.logger.log.info("[Biomarker outlier detection] - [Sample location: '{}'] - Using following outlier detection methods: {}".format(
            sample_location, self.biomarker_outlier_statistics))
        stat_dict = dict()
        progress_bar = self.logger.get_progress_bar(CalculatedColumns.get_num_of_unprocessed(measurements_df), "Analyzing biomarker outliers")
        for index, current_measurement in measurements_df.iterrows():
            if CalculatedColumns.needs_processing(current_measurement):
                progress_bar.update(1)
                for biomarker1, biomarker2 in itertools.combinations(Columns.get_biomarker_columns(), 2):
                    add_default_biomarker_statistic(stat_dict, biomarker1, biomarker2)
                    biomarker_ratio = biomarker1 + "/" + biomarker2
                    if current_measurement[biomarker_ratio] and not math.isnan(current_measurement[biomarker_ratio]):  # only if biomarker ratio is not None for current measurement
                        stat_dict[biomarker_ratio]["total"] += 1
                        last_biomarker_ratios = self.__get_previous_biomarkers_ratios(measurements_df, index, biomarker1, biomarker2)
                        if len(last_biomarker_ratios) < self.min_number_biomarkers_for_outlier_detection:
                            SewageFlag.add_flag_to_index_column(measurements_df, index, CalculatedColumns.get_biomaker_ratio_flag(biomarker1, biomarker2),
                                                                SewageFlag.NOT_ENOUGH_PREVIOUS_BIOMARKER_VALUES)
                            self.logger.log.debug(
                                "[Biomarker outlier detection] - [Sample location: '{}'] - [Collection date: '{}'] - "
                                "Skipping biomarker outlier detection. Found '{}' previous measurements. "
                                "Minimal '{}' previous measurements are required.".format(sample_location, current_measurement[Columns.DATE.value],
                                                                                          len(last_biomarker_ratios),
                                                                                          self.min_number_biomarkers_for_outlier_detection))
                            stat_dict[biomarker1 + "/" + biomarker2]["skipped"] += 1
                            continue
                        is_outlier = detect_outliers(self.biomarker_outlier_statistics, last_biomarker_ratios[biomarker_ratio], current_measurement[biomarker_ratio])
                        if is_outlier:
                            stat_dict[biomarker_ratio]["failed"] += 1
                            SewageFlag.add_flag_to_index_column(measurements_df, index, CalculatedColumns.get_biomaker_ratio_flag(biomarker1, biomarker2), SewageFlag.BIOMARKER_RATIO_OUTLIER)
                        else:
                            stat_dict[biomarker1 + "/" + biomarker2]["passed"] += 1
                    else:
                        stat_dict[biomarker1 + "/" + biomarker2]["skipped"] += 1
                        self.logger.log.debug(
                            "[Biomarker outlier detection] - [Sample location: '{}'] - [Collection date: '{}'] - "
                            "Empty biomarker ratio '{}' with value: '{}'. Biomarker ratio will be skipped."
                            .format(sample_location, current_measurement[Columns.DATE.value],
                                    biomarker_ratio,
                                    current_measurement[biomarker_ratio]))
        progress_bar.close()
        if self.plotting:
            plot_biomarker_outlier_summary(measurements_df, sample_location, os.path.join(self.output_folder, "plots", "biomarker"),
                                           self.biomarker_outlier_statistics, self.interactive)
        self.logger.log.info("[Biomarker outlier detection] - [Sample location: '{}'] - \n"
                             "Biomarkers ratios:\n{}".format(sample_location,
                                                             pretty_print_biomarker_statistic(stat_dict)))

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

    def __flag_outliers_with_all_ratio_outliers(self, biomarker_outlier_set: set, index: int, measurements_df: pd.DataFrame) -> None:
        for biomarker in biomarker_outlier_set:
            SewageFlag.add_flag_to_index_column(measurements_df, index, CalculatedColumns.get_biomarker_flag(biomarker),
                                                SewageFlag.BIOMARKER_VALIDATED_OUTLIER)

    def assign_biomarker_outliers_based_on_ratio_flags(self, sample_location, measurements_df: pd.DataFrame):
        for index, current_measurement in measurements_df.iterrows():
            if CalculatedColumns.needs_processing(current_measurement):
                usable_biomarkers, biomarker_ratio_flags = self.__get_usable_biomarkers(current_measurement)
                all_ratios_outlier_biomarker_set = self.__mark_biomakers_with_all_outliers(current_measurement)
                # if all ratios for a given biomarker are flagged as outliers then the biomarkers is flagged as validated outlier
                self.__flag_outliers_with_all_ratio_outliers(all_ratios_outlier_biomarker_set, index, measurements_df)
                for biomarker1, biomarker2 in itertools.combinations(usable_biomarkers, 2):
                    # is biomarker1 or biomarker2 is not already flagged as outlier before
                    if not biomarker1 in all_ratios_outlier_biomarker_set and not biomarker2 in all_ratios_outlier_biomarker_set:
                        is_outlier = SewageFlag.is_flag(current_measurement[CalculatedColumns.get_biomaker_ratio_flag(biomarker1, biomarker2)], SewageFlag.BIOMARKER_RATIO_OUTLIER)
                        if is_outlier:
                            SewageFlag.add_flag_to_index_column(measurements_df, index, CalculatedColumns.get_biomarker_flag(biomarker1), SewageFlag.BIOMARKER_PROBABLE_OUTLIER)
                            SewageFlag.add_flag_to_index_column(measurements_df, index, CalculatedColumns.get_biomarker_flag(biomarker2), SewageFlag.BIOMARKER_PROBABLE_OUTLIER)

    def analyze_usable_biomarkers(self, sample_location, measurements_df: pd.DataFrame):
        """
        Sums the number of biomarkers which are not marked as outliers using the outlier columns.
        Flags measurments which have less than specified number of valid biomarkers.
        """
        for index, current_measurement in measurements_df.iterrows():
            if CalculatedColumns.needs_processing(current_measurement):
                num_usable_biomarkers = 0
                is_probable_outlier_flag = False
                for biomarker in Columns.get_biomarker_columns():
                    is_below_threshold_or_empty = SewageFlag.is_flag(current_measurement[CalculatedColumns.get_biomarker_flag(biomarker)], SewageFlag.BIOMARKER_BELOW_THRESHOLD_OR_EMPTY)
                    is_validated_outlier = SewageFlag.is_flag(current_measurement[CalculatedColumns.get_biomarker_flag(biomarker)], SewageFlag.BIOMARKER_VALIDATED_OUTLIER)
                    is_probable_outlier = SewageFlag.is_flag(current_measurement[CalculatedColumns.get_biomarker_flag(biomarker)], SewageFlag.BIOMARKER_PROBABLE_OUTLIER)
                    if not is_below_threshold_or_empty and not is_validated_outlier:
                        num_usable_biomarkers += 1
                    if is_probable_outlier:
                        is_probable_outlier_flag = True
                measurements_df.at[index, CalculatedColumns.NUMBER_OF_USABLE_BIOMARKERS.value] = num_usable_biomarkers
                #            if num_usable_biomarkers < 2:  # check if less than 2 biomarkers are available
                #                SewageFlag.add_flag_to_index_column(measurements_df, index, Columns.FLAG.value, SewageFlag.MIN_BIOMARKER_NUMBER_NOT_REACHED)
                if is_probable_outlier_flag:  # in case any biomarker was marked as probable outlier set flag for normalization
                    SewageFlag.add_flag_to_index_column(measurements_df, index, CalculatedColumns.FLAG.value, SewageFlag.BIOMARKER_PROBABLE_OUTLIER)
        if CalculatedColumns.NUMBER_OF_USABLE_BIOMARKERS.value in measurements_df:
            biomarker_stat = measurements_df[CalculatedColumns.NUMBER_OF_USABLE_BIOMARKERS.value].value_counts().sort_index(ascending=True)
            stat_output = ""
            for num_biomarker, num_samples in zip(biomarker_stat.index.values, biomarker_stat.values):
                stat_output += "\t{} usable biomarkers found in {} samples\n".format(int(num_biomarker), num_samples)
            self.logger.log.info("[Select Biomarkers] - [Sample location: '{}'] - \n"
                                 "{}".format(sample_location, stat_output))

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
