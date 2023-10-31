# Created by alex at 03.07.23
import pandas as pd
import math
from .utils import *
from .statistics import *
from .plotting import *


class SewageNormalization:
    def __init__(self, sewageStat: SewageStat, max_number_of_flags_for_outlier: int,
                 min_number_of_biomarkers_for_normalization: int,
                 base_reproduction_value_factor: float,
                 num_previous_days_reproduction_factor: int,
                 output_folder: str):
        self.sewageStat = sewageStat
        self.max_number_of_flags_for_outlier = max_number_of_flags_for_outlier
        self.min_number_of_biomarkers_for_normalization = min_number_of_biomarkers_for_normalization
        self.base_reproduction_value_factor = base_reproduction_value_factor
        self.num_previous_days_reproduction_factor = num_previous_days_reproduction_factor
        self.output_folder = output_folder
        self.logger = SewageLogger(self.output_folder)

    def __get_usable_biomarkers(self, current_measurement: pd.Series) -> []:
        biomarker_values = []
        for biomarker in Columns.get_biomarker_columns():
            if SewageFlag.is_not_flag(current_measurement[CalculatedColumns.get_biomarker_flag(biomarker)], SewageFlag.BIOMARKER_VALIDATED_OUTLIER) and \
                    SewageFlag.is_not_flag(current_measurement[CalculatedColumns.get_biomarker_flag(biomarker)], SewageFlag.BIOMARKER_BELOW_THRESHOLD_OR_EMPTY):
                biomarker_value = current_measurement[biomarker]
                if biomarker_value and not math.isnan(biomarker_value):
                    biomarker_values.append(biomarker_value)
        return biomarker_values

    def __normalize_with_sewage_flow(self, measurements_df: pd.DataFrame, index):
        current_measurement = measurements_df.iloc[index]
        if current_measurement[CalculatedColumns.NUMBER_OF_USABLE_BIOMARKERS.value] >= self.min_number_of_biomarkers_for_normalization:
            if SewageFlag.is_not_flag(current_measurement[CalculatedColumns.FLAG.value], SewageFlag.MISSING_MEAN_SEWAGE_FLOW):
                biomarker_values = self.__get_usable_biomarkers(current_measurement)
                mean_sewage_flow = (current_measurement[Columns.MEAN_SEWAGE_FLOW.value] / 1000) * 60 * 60 * 24  # from l/s --> m³/day
                mean_biomarker_value = np.mean(biomarker_values) * 1000 * 1000  # from genecopies/ml --> mean genecopies/m³
                normalized_mean_biomarker = round(mean_biomarker_value * mean_sewage_flow, 2)
                return normalized_mean_biomarker
        else:
            SewageFlag.add_flag_to_index_column(measurements_df, index, CalculatedColumns.FLAG.value, SewageFlag.NOT_ENOUGH_BIOMARKERS_FOR_NORMALIZATION)
        return None

    def __detect_basic_reproduction_number_outliers(self, sample_location: str, measurements_df: pd.DataFrame, index) -> None:
        current_measurement = measurements_df.iloc[index]
        last_values_one_week = get_last_N_month_and_days(measurements_df, current_measurement, CalculatedColumns.NORMALIZED_MEAN_BIOMARKERS.value,
                                                         num_month=0, num_days=self.num_previous_days_reproduction_factor,
                                                         sewage_flag=SewageFlag.REPRODUCTION_NUMBER_OUTLIER)
        last_values_one_week = last_values_one_week[last_values_one_week[CalculatedColumns.NORMALIZED_MEAN_BIOMARKERS.value] > 0]
        if last_values_one_week.shape[0] > 0:
            current_mean_normalized_biomarker = current_measurement[CalculatedColumns.NORMALIZED_MEAN_BIOMARKERS.value]
            last_mean_normalized_biomarker = np.mean(last_values_one_week[CalculatedColumns.NORMALIZED_MEAN_BIOMARKERS.value])
            if last_mean_normalized_biomarker > 0 and current_mean_normalized_biomarker > 0:  # no division by zero
                reproduction_factor = current_mean_normalized_biomarker / last_mean_normalized_biomarker
                measurements_df.at[index, CalculatedColumns.BASE_REPRODUCTION_FACTOR.value] = reproduction_factor
                if reproduction_factor > self.base_reproduction_value_factor or reproduction_factor < (1 / self.base_reproduction_value_factor):
                    SewageFlag.add_flag_to_index_column(measurements_df, index, CalculatedColumns.FLAG.value, SewageFlag.REPRODUCTION_NUMBER_OUTLIER)
                    self.sewageStat.add_reproduction_factor_outlier('failed')
            else:
                SewageFlag.add_flag_to_index_column(measurements_df, index, CalculatedColumns.FLAG.value, SewageFlag.REPRODUCTION_NUMBER_OUTLIER_SKIPPED)
                self.sewageStat.add_reproduction_factor_outlier('skipped')
        else:
            SewageFlag.add_flag_to_index_column(measurements_df, index, CalculatedColumns.FLAG.value, SewageFlag.REPRODUCTION_NUMBER_OUTLIER_SKIPPED)
            self.sewageStat.add_reproduction_factor_outlier('skipped')

    def normalize_biomarker_values(self, sample_location: str, measurements: pd.DataFrame, index) -> None:
        normalized_mean_biomarker = self.__normalize_with_sewage_flow(measurements, index)
        if normalized_mean_biomarker:
            measurements.at[index, CalculatedColumns.NORMALIZED_MEAN_BIOMARKERS.value] = normalized_mean_biomarker
            self.__detect_basic_reproduction_number_outliers(sample_location, measurements, index)
        else:
            # no normalized mean biomarker could be calculated --> too less usable biomarkers
            self.sewageStat.add_reproduction_factor_outlier('failed')
            SewageFlag.add_flag_to_index_column(measurements, index, CalculatedColumns.FLAG.value, SewageFlag.REPRODUCTION_NUMBER_OUTLIER_SKIPPED)


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
            is_outlier = SewageFlag.is_flag(current_measurement[CalculatedColumns.FLAG.value], CalculatedColumns.get_surrogate_outlier_flag(sVirus))
            if is_outlier:
                num_surrogate_virus_flags += 1
            all_surrogate_virus_outliers.append(is_outlier)
        return all(all_surrogate_virus_outliers), num_surrogate_virus_flags

    def __is_sewage_flow_outlier(self, current_measurement: pd.Series) -> bool:
        flag_value = current_measurement[CalculatedColumns.FLAG.value]
        # too heavy precipitation or typo in sewage flow
        return SewageFlag.is_flag(flag_value, SewageFlag.SEWAGE_FLOW_PRECIPITATION)

    def __is_sewage_heavy_precipitation_outlier(self, current_measurement: pd.Series) -> bool:
        flag_value = current_measurement[CalculatedColumns.FLAG.value]
        # too heavy precipitation or typo in sewage flow
        return SewageFlag.is_flag(flag_value, SewageFlag.SEWAGE_FLOW_HEAVY_PRECIPITATION)

    def __is_sewage_probable_typo(self, current_measurement: pd.Series) -> bool:
        flag_value = current_measurement[CalculatedColumns.FLAG.value]
        return SewageFlag.is_flag(flag_value, SewageFlag.SEWAGE_FLOW_PROBABLE_TYPO)

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


    def __clear_biomarker_ratio_outliers(self, measurements: pd.DataFrame, index):
        for biomarker1, biomarker2 in itertools.combinations(Columns.get_biomarker_columns(), 2):
            current_flag = measurements.iloc[index][CalculatedColumns.get_biomaker_ratio_flag(biomarker1, biomarker2)]
            if SewageFlag.is_flag(current_flag, SewageFlag.BIOMARKER_RATIO_OUTLIER):
                SewageFlag.add_flag_to_index_column(measurements, index, CalculatedColumns.get_biomaker_ratio_flag(biomarker1, biomarker2), SewageFlag.BIOMARKER_RATIO_OUTLIER_REMOVED)


    def decide_biomarker_usable_based_on_flags(self, sample_location: str, measurements: pd.DataFrame, index):
        current_measurement = measurements.iloc[index]
        usable = 0
        num_flags = 0
        is_outlier = False
        if self.__are_comments_not_empty(current_measurement):
            num_flags += 1
        #  less than 2 biomarker values available
        if self.__are_too_less_biomarkers(current_measurement):
            is_outlier = True
            self.sewageStat.add_outliers('Min num biomarkers not reached')
            CalculatedColumns.add_outlier_reason(measurements, index, 'Min num biomarkers not reached')
        num_flags += self.__num_biomarkers_flagged(current_measurement)
        # surrogate_virus_flags
        are_both_surrogate_virus_outliers, num_surrogate_virus_flags = self.__are_both_surrogate_virus_outliers(current_measurement)
        if are_both_surrogate_virus_outliers:
            self.sewageStat.add_outliers('Surrogate virus outlier')
            CalculatedColumns.add_outlier_reason(measurements, index, 'Surrogate virus outlier')
            is_outlier = True
        else:
            num_flags += num_surrogate_virus_flags
        if self.__is_sewage_flow_outlier(current_measurement):   # precipitation outlier
            if num_surrogate_virus_flags > 0:   #
                self.sewageStat.add_outliers('Sewage flow outlier')
                CalculatedColumns.add_outlier_reason(measurements, index, 'Sewage flow outlier')
                is_outlier = True
        if self.__is_sewage_heavy_precipitation_outlier(current_measurement):   # heavy precipitation outlier
            self.sewageStat.add_outliers('Sewage flow outlier')
            CalculatedColumns.add_outlier_reason(measurements, index, 'Sewage flow outlier')
            is_outlier = True
        if self.__is_sewage_probable_typo(current_measurement):   # probable typo outlier
            self.sewageStat.add_outliers('Sewage flow outlier')
            CalculatedColumns.add_outlier_reason(measurements, index, 'Sewage flow outlier')
            is_outlier = True
        num_flags += self.__get_num_water_quality_flags(current_measurement)
        if self.__is_reproduction_factor_outlier(current_measurement):
            self.sewageStat.add_outliers('Reproduction factor outlier')
            CalculatedColumns.add_outlier_reason(measurements, index, 'Reproduction factor outlier')
            is_outlier = True
        if num_flags > self.max_number_of_flags_for_outlier:
            CalculatedColumns.add_outlier_reason(measurements, index, 'Too many flags')
            self.sewageStat.add_outliers('Too many flags')
            is_outlier = True
        if is_outlier:
            measurements.at[index, CalculatedColumns.USABLE.value] = False
        else:
            usable += 1
            measurements.at[index, CalculatedColumns.USABLE.value] = True
            # remove biomarker ratio outliers in case the sample is valid
            self.__clear_biomarker_ratio_outliers(measurements, index)


