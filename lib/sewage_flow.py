# Created by alex at 03.07.23
import math
from .plotting import *
from .statistics import *
from .utils import *


class SewageFlow:

    def __init__(self, output_folder, sewageStat: SewageStat, sewage_plants2dry_weather_flow: dict, fraction_last_samples_for_dry_flow: float,
                 min_num_samples_for_mean_dry_flow: int, heavy_precipitation_factor: int, mean_sewage_flow_below_typo_factor: float,
                 mean_sewage_flow_above_typo_factor: float):
        self.output_folder = output_folder
        self.sewageStat = sewageStat
        self.sewage_plants2dry_weather_flow = sewage_plants2dry_weather_flow
        self.fraction_last_samples_for_dry_flow = fraction_last_samples_for_dry_flow
        self.min_num_samples_for_mean_dry_flow = min_num_samples_for_mean_dry_flow
        self.heavy_precipitation_factor = heavy_precipitation_factor
        self.mean_sewage_flow_below_typo_factor = mean_sewage_flow_below_typo_factor
        self.mean_sewage_flow_above_typo_factor = mean_sewage_flow_above_typo_factor
        self.logger = SewageLogger(self.output_folder)

    def sewage_flow_quality_control(self, sample_location, measurements: pd.DataFrame, index):
        dry_flow, is_mean_dry_weather_flow_available = self.__is_mean_flow_above_dry_flow(sample_location, measurements, index)

    def __get_mean_flow_based_on_last_min_values(self, sample_location, measurements_df: pd.DataFrame, current_measurement):
        last_measurements = measurements_df[measurements_df[Columns.DATE.value] < current_measurement[Columns.DATE.value]]
        # omit samples where any sewage flow tag was set
        last_measurements = last_measurements[SewageFlag.is_not_flag_set_for_series(last_measurements[CalculatedColumns.FLAG.value], SewageFlag.SEWAGE_FLOW_HEAVY_PRECIPITATION)]
        last_measurements = last_measurements[SewageFlag.is_not_flag_set_for_series(last_measurements[CalculatedColumns.FLAG.value], SewageFlag.SEWAGE_FLOW_PROBABLE_TYPO)]
        last_measurements = last_measurements[SewageFlag.is_not_flag_set_for_series(last_measurements[CalculatedColumns.FLAG.value], SewageFlag.MISSING_MEAN_SEWAGE_FLOW)]
        mean_sewage_flows = last_measurements[Columns.MEAN_SEWAGE_FLOW.value]
        if mean_sewage_flows.shape[0] < self.min_num_samples_for_mean_dry_flow:
            self.logger.log.debug("[Sewage flow] - [Sample location: '{}'] - Less than '{}' "
                                  "previous samples obtained. Skipping sewage flow QC...".format(sample_location, self.min_num_samples_for_mean_dry_flow))
            return None
        # round up to next integer
        num_last_N_samples = math.ceil(mean_sewage_flows.shape[0] * self.fraction_last_samples_for_dry_flow)
        smallest_sewage_flows = mean_sewage_flows.nsmallest(num_last_N_samples, keep='all')
        mean_dry_flow_estimation = np.mean(smallest_sewage_flows)
        return mean_dry_flow_estimation

    def __get_dry_flow(self, sample_location, measurements_df: pd.DataFrame, current_measurement):
        is_mean_dry_weather_flow_available = False
        mean_dry_flow_estimation = self.__get_mean_flow_based_on_last_min_values(sample_location, measurements_df, current_measurement)
        if mean_dry_flow_estimation:
            is_mean_dry_weather_flow_available = True
            return mean_dry_flow_estimation, is_mean_dry_weather_flow_available
        else:
            if sample_location in self.sewage_plants2dry_weather_flow:
                dry_flow = self.sewage_plants2dry_weather_flow[sample_location]
                if dry_flow:
                    return dry_flow, is_mean_dry_weather_flow_available
        return None, is_mean_dry_weather_flow_available

    def __is_mean_flow_above_dry_flow(self, sample_location, measurements: pd.DataFrame, index):
        is_mean_dry_weather_flow_available, dry_flow = False, None
        current_measurement = measurements.iloc[index]
        if SewageFlag.is_not_flag(current_measurement[CalculatedColumns.FLAG.value], SewageFlag.MISSING_MEAN_SEWAGE_FLOW):
            dry_flow, is_mean_dry_weather_flow_available = self.__get_dry_flow(sample_location, measurements, current_measurement)
            if dry_flow:
                current_mean_flow = current_measurement[Columns.MEAN_SEWAGE_FLOW.value]
                # is the mean flow a factor of N (default: 9) higher than the dry weather flow --> probable typo
                if (current_mean_flow / dry_flow) > self.mean_sewage_flow_above_typo_factor:
                    self.sewageStat.add_sewage_flow_outlier('probable_typo')
                    measurements.at[index, CalculatedColumns.FLAG.value] += SewageFlag.SEWAGE_FLOW_PROBABLE_TYPO.value
                # is the mean flow a factor of N (default: 2) higher than the dry weather flow --> heavy precipitation
                elif (current_mean_flow / dry_flow) > self.heavy_precipitation_factor:
                    self.sewageStat.add_sewage_flow_outlier('heavy_precipitation')
                    measurements.at[index, CalculatedColumns.FLAG.value] += SewageFlag.SEWAGE_FLOW_HEAVY_PRECIPITATION.value
                # is the mean flow a factor of N (default: 1.5) less than the dry weather flow --> probably typo
                elif (dry_flow / current_mean_flow) > self.mean_sewage_flow_below_typo_factor:
                    self.sewageStat.add_sewage_flow_outlier('probable_typo')
                    measurements.at[index, CalculatedColumns.FLAG.value] += SewageFlag.SEWAGE_FLOW_PROBABLE_TYPO.value
                else:
                    self.sewageStat.add_sewage_flow_outlier('passed')
            else:
                measurements.at[index, CalculatedColumns.FLAG.value] += SewageFlag.SEWAGE_FLOW_NOT_ENOUGH_PREVIOUS_VALUES.value
                self.sewageStat.add_sewage_flow_outlier('skipped')
        else:
            measurements.at[index, CalculatedColumns.FLAG.value] += SewageFlag.SEWAGE_FLOW_NOT_ENOUGH_PREVIOUS_VALUES.value
            self.sewageStat.add_sewage_flow_outlier('skipped')
        return dry_flow, is_mean_dry_weather_flow_available

