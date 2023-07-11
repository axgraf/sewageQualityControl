# Created by alex at 03.07.23
import pandas as pd
import math
from .utils import *
from .plotting import *


class SewageNormalization:
    def __init__(self, max_number_of_flags_for_outlier: int,
                 min_number_of_biomarkers_for_normalization: int,
                 base_reproduction_value_factor: float,
                 output_folder: str):
        self.max_number_of_flags_for_outlier = max_number_of_flags_for_outlier
        self.min_number_of_biomarkers_for_normalization = min_number_of_biomarkers_for_normalization
        self.base_reproduction_value_factor = base_reproduction_value_factor
        self.output_folder = output_folder
        self.logger = SewageLogger(self.output_folder)

    def __get_usable_biomarkers(self, current_measurement: pd.Series) -> []:
        biomarker_values = []
        for biomarker in Columns.get_biomarker_columns():
            if SewageFlag.is_not_flag(current_measurement[CalculatedColumns.get_biomarker_flag(biomarker)], SewageFlag.BIOMARKER_VALIDATED_OUTLIER) or \
                    SewageFlag.is_not_flag(current_measurement[CalculatedColumns.get_biomarker_flag(biomarker)], SewageFlag.BIOMARKER_BELOW_THRESHOLD_OR_EMPTY):
                biomarker_value = current_measurement[biomarker]
                if biomarker_value and not math.isnan(biomarker_value):
                    biomarker_values.append(biomarker_value)
        return biomarker_values

    def __normalize_with_sewage_flow(self, measurements_df: pd.DataFrame, index: int) -> float:
        current_measurement = measurements_df.iloc[index]
        if current_measurement[CalculatedColumns.NUMBER_OF_USABLE_BIOMARKERS.value] >= self.min_number_of_biomarkers_for_normalization:
            if SewageFlag.is_not_flag(current_measurement[CalculatedColumns.FLAG.value], SewageFlag.MISSING_MEAN_SEWAGE_FLOW):
                biomarker_values = self.__get_usable_biomarkers(current_measurement)
                if len(biomarker_values) == 0:
                    print("here")
                mean_sewage_flow = current_measurement[Columns.MEAN_SEWAGE_FLOW.value]
                mean_biomarker_value = np.mean(biomarker_values)
                normalized_mean_biomarker = round(mean_biomarker_value * mean_sewage_flow, 2)
                return normalized_mean_biomarker
        else:
            SewageFlag.add_flag_to_index_column(measurements_df, index, CalculatedColumns.FLAG.value, SewageFlag.NOT_ENOUGH_BIOMARKERS_FOR_NORMALIZATION)
        return None

    def __detect_basic_reproduction_number_outliers(self, sample_location: str, measurements_df: pd.DataFrame, index: int) -> None:
        current_measurement = measurements_df.iloc[index]
        num_previous_days = 7
        last_values_one_week = get_last_N_month_and_days(measurements_df, current_measurement, CalculatedColumns.NORMALIZED_MEAN_BIOMARKERS.value, num_month=0, num_days=num_previous_days,
                                                         sewage_flag=None)
        if last_values_one_week.shape[0] > 0:
            last_mean_normalized_biomarker = np.mean(last_values_one_week[CalculatedColumns.NORMALIZED_MEAN_BIOMARKERS.value])
            current_mean_normalized_biomarker = current_measurement[CalculatedColumns.NORMALIZED_MEAN_BIOMARKERS.value]
            if last_mean_normalized_biomarker > 0 and current_mean_normalized_biomarker > 0:  # no division by zero
                reproduction_factor = current_mean_normalized_biomarker / last_mean_normalized_biomarker
                measurements_df.at[index, CalculatedColumns.BASE_REPRODUCTION_FACTOR.value] = reproduction_factor
                if reproduction_factor > self.base_reproduction_value_factor or reproduction_factor < (1 / self.base_reproduction_value_factor):
                    SewageFlag.add_flag_to_index_column(measurements_df, index, CalculatedColumns.FLAG.value, SewageFlag.REPRODUCTION_NUMBER_OUTLIER)
            else:
                SewageFlag.add_flag_to_index_column(measurements_df, index, CalculatedColumns.FLAG.value, SewageFlag.REPRODUCTION_NUMBER_OUTLIER_SKIPPED)
        else:
            SewageFlag.add_flag_to_index_column(measurements_df, index, CalculatedColumns.FLAG.value, SewageFlag.REPRODUCTION_NUMBER_OUTLIER_SKIPPED)
            self.logger.log.debug("[Biomarker reproduction number outlier] - [Sample location: '{}'] - [Collection date: '{}'] - "
                                  "No previous mean normalized biomarker values found within the last {} days. "
                                  "Skipping reproduction outlier analysis... ".format(sample_location, current_measurement[Columns.DATE.value], num_previous_days))

    def normalize_biomarker_values(self, sample_location: str, measurements_df: pd.DataFrame) -> None:
        for index, current_measurement in measurements_df.iterrows():
            if CalculatedColumns.needs_processing(current_measurement):
                if index == 21:
                    print("here")
                normalized_mean_biomarker = self.__normalize_with_sewage_flow(measurements_df, index)
                if normalized_mean_biomarker:
                    measurements_df.at[index, CalculatedColumns.NORMALIZED_MEAN_BIOMARKERS.value] = normalized_mean_biomarker
                    self.__detect_basic_reproduction_number_outliers(sample_location, measurements_df, index)
                else:
                    ## no normalized mean biomarker could be calculated --> too less usable biomarkers
                    SewageFlag.add_flag_to_index_column(measurements_df, index, CalculatedColumns.FLAG.value, SewageFlag.REPRODUCTION_NUMBER_OUTLIER_SKIPPED)
                    self.logger.log.warn("[Biomarker normalization] - [Sample location: '{}'] -  "
                                         "Normalized biomarker could not be calculated. Either too less valid biomarkers or mean sewage flow is missing".
                                         format(sample_location))
        plot_biomarker_normalization(measurements_df, sample_location, os.path.join(self.output_folder, "plots", "biomarker_normalization"))


    def __are_comments_not_empty(self, current_measurement: pd.Series) -> bool:
        flag_value = current_measurement[CalculatedColumns.FLAG.value]
        return SewageFlag.is_flag(flag_value, SewageFlag.COMMENT_NOT_EMPTY)

    def __are_too_less_biomarkers(self, current_measurement: pd.Series) -> bool:
        flag_value = current_measurement[CalculatedColumns.FLAG.value]
        return SewageFlag.is_flag(flag_value, SewageFlag.NOT_ENOUGH_BIOMARKERS_FOR_NORMALIZATION)

    def __num_biomarkers_flagged(self, current_measurement: pd.Series) -> int:
        num_biomarker_flags = 0
        for biomarker in Columns.get_biomarker_columns():
            flag_value = current_measurement[CalculatedColumns.get_biomarker_flag(biomarker)]
            if SewageFlag.is_flag(flag_value, SewageFlag.BIOMARKER_PROBABLE_OUTLIER) or \
                    SewageFlag.is_flag(flag_value, SewageFlag.BIOMARKER_VALIDATED_OUTLIER):
                num_biomarker_flags += 1
        return num_biomarker_flags

    def __are_both_surrogate_virus_outliers(self, current_measurement: pd.Series) -> (bool, int):
        all_surrogate_virus_outliers = []
        num_surrogate_virus_flags = 0
        for sVirus in Columns.get_surrogatevirus_columns():
            is_outlier = SewageFlag.is_flag(current_measurement[CalculatedColumns.get_surrogate_outlier_flag(sVirus)], SewageFlag.SURROGATEVIRUS_OUTLIER)
            if is_outlier:
                num_surrogate_virus_flags += 1
            all_surrogate_virus_outliers.append(is_outlier)
        return all(all_surrogate_virus_outliers), num_surrogate_virus_flags

    def __is_sewage_flow_outlier(self, current_measurement: pd.Series) -> bool:
        flag_value = current_measurement[CalculatedColumns.FLAG.value]
        # too heavy precipitation or typo in sewage flow
        return SewageFlag.is_flag(flag_value, SewageFlag.SEWAGE_FLOW_PROBABLE_TYPO) or \
            SewageFlag.is_flag(flag_value, SewageFlag.SEWAGE_FLOW_HEAVY_PRECIPITATION)

    def __get_num_water_quality_flags(self, current_measurement: pd.Series) -> int:
        flag_value = current_measurement[CalculatedColumns.FLAG.value]
        num_water_quality_flags = 0
        if SewageFlag.is_flag(flag_value, SewageFlag.AMMONIUM_OUTLIER):
            num_water_quality_flags += 1
        if SewageFlag.is_flag(flag_value, SewageFlag.CONDUCTIVITY_OUTLIER):
            num_water_quality_flags += 1
        return num_water_quality_flags

    def __is_reproduction_factor_outlier(self, current_measurement: pd.Series) -> bool:
        flag_value = current_measurement[CalculatedColumns.FLAG.value]
        return SewageFlag.is_flag(flag_value, SewageFlag.REPRODUCTION_NUMBER_OUTLIER)

    def decide_biomarker_usable_based_on_flags(self, sample_location: str, measurements_df: pd.DataFrame):
        usable = 0
        outlier_stat_dict = dict()
        outlier_stat_dict.setdefault('Min num biomarkers not reached', 0)
        outlier_stat_dict.setdefault('Surrogate virus outlier', 0)
        outlier_stat_dict.setdefault('Sewage flow outlier', 0)
        outlier_stat_dict.setdefault('Reproduction factor outlier', 0)
        outlier_stat_dict.setdefault('Too many flags', 0)

        for index, current_measurement in measurements_df.iterrows():
            if CalculatedColumns.needs_processing(current_measurement):
                num_flags = 0
                is_outlier = False
                if self.__are_comments_not_empty(current_measurement):
                    num_flags += 1
                #  less than 2 biomarker values available
                if self.__are_too_less_biomarkers(current_measurement):
                    is_outlier = True
                    outlier_stat_dict['Min num biomarkers not reached'] += 1
                    CalculatedColumns.add_outlier_reason(measurements_df, index, 'Min num biomarkers not reached')
                num_flags += self.__num_biomarkers_flagged(current_measurement)
                # surrogate_virus_flags
                are_both_surrogate_virus_outliers, num_surrogate_virus_flags = self.__are_both_surrogate_virus_outliers(current_measurement)
                num_flags += num_surrogate_virus_flags
                if are_both_surrogate_virus_outliers:
                    outlier_stat_dict['Surrogate virus outlier'] += 1
                    CalculatedColumns.add_outlier_reason(measurements_df, index, 'Surrogate virus outlier')
                    is_outlier = True
                if self.__is_sewage_flow_outlier(current_measurement):
                    outlier_stat_dict['Sewage flow outlier'] += 1
                    CalculatedColumns.add_outlier_reason(measurements_df, index, 'Sewage flow outlier')
                    is_outlier = True
                num_flags += self.__get_num_water_quality_flags(current_measurement)
                if self.__is_reproduction_factor_outlier(current_measurement):
                    outlier_stat_dict['Reproduction factor outlier'] += 1
                    CalculatedColumns.add_outlier_reason(measurements_df, index, 'Reproduction factor outlier')
                    is_outlier = True
                if num_flags > self.max_number_of_flags_for_outlier:
                    CalculatedColumns.add_outlier_reason(measurements_df, index, 'Too many flags')
                    outlier_stat_dict['Too many flags'] += 1
                    is_outlier = True
                if is_outlier:
                    measurements_df.at[index, CalculatedColumns.USABLE.value] = False
                else:
                    usable += 1
                    measurements_df.at[index, CalculatedColumns.USABLE.value] = True
        num2process = CalculatedColumns.get_num_of_unprocessed(measurements_df)
        self.logger.log.info("[Determine outliers based on all flags] - [Sample location: '{}'] -  "
                             "Detected outliers '{}/{}' from all quality control steps.\nOutlier reasons:\n{}".
                             format(sample_location, (num2process - usable), num2process, outlier_stat_dict))
        plot_general_outliers(measurements_df, sample_location,  os.path.join(self.output_folder, "plots", "outlier"))
        print("here")


