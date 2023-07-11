# Created by alex at 03.07.23
import math
from .plotting import *
from .utils import *
from .database import SewageDatabase


class SewageFlow:

    def __init__(self, output_folder, interactive, sewage_plants2dry_weather_flow: dict, fraction_last_samples_for_dry_flow: float,
                 min_num_samples_for_mean_dry_flow: int, heavy_precipitation_factor: int, mean_sewage_flow_below_typo_factor: float,
                 mean_sewage_flow_above_typo_factor: float):
        self.output_folder = output_folder
        self.interactive = interactive
        self.sewage_plants2dry_weather_flow = sewage_plants2dry_weather_flow
        self.fraction_last_samples_for_dry_flow = fraction_last_samples_for_dry_flow
        self.min_num_samples_for_mean_dry_flow = min_num_samples_for_mean_dry_flow
        self.heavy_precipitation_factor = heavy_precipitation_factor
        self.mean_sewage_flow_below_typo_factor = mean_sewage_flow_below_typo_factor
        self.mean_sewage_flow_above_typo_factor = mean_sewage_flow_above_typo_factor
        self.logger = SewageLogger(self.output_folder)
        self.database = SewageDatabase(self.output_folder)

    def sewage_flow_quality_control(self, sample_location, measurements_df: pd.DataFrame):
        dry_flow, is_mean_dry_weather_flow_available = self.__is_mean_flow_above_dry_flow(sample_location, measurements_df)
        plot_sewage_flow(measurements_df, sample_location, os.path.join(self.output_folder, "plots", "sewage_flow"),
                         is_mean_dry_weather_flow_available, dry_flow, self.interactive)

    def __get_mean_flow_based_on_last_min_values(self, sample_location, measurements_df: pd.DataFrame, current_measurement):
        last_measurements = measurements_df[measurements_df[Columns.DATE.value] < current_measurement[Columns.DATE.value]]
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

    def __is_mean_flow_above_dry_flow(self, sample_location, measurements_df: pd.DataFrame):
        is_mean_dry_weather_flow_available, dry_flow = False, None
        mean_flow_quality_stats = dict()
        for k in ['passed', 'skipped', 'probable_typo', 'heavy_precipitation', 'total']:
            mean_flow_quality_stats.setdefault(k, 0)
        progress_bar = self.logger.get_progress_bar(CalculatedColumns.get_num_of_unprocessed(measurements_df), "Analyzing mean sewage flows")
        for index, current_measurement in measurements_df.iterrows():
            if CalculatedColumns.needs_processing(current_measurement):
                progress_bar.update(1)
                dry_flow, is_mean_dry_weather_flow_available = self.__get_dry_flow(sample_location, measurements_df, current_measurement)
                if dry_flow:
                    mean_flow_quality_stats['total'] += 1
                    current_mean_flow = current_measurement[Columns.MEAN_SEWAGE_FLOW.value]
                    # is the mean flow a factor of N (default: 9) higher than the dry weather flow --> probable typo
                    if (current_mean_flow / dry_flow) > self.mean_sewage_flow_above_typo_factor:
                        mean_flow_quality_stats['probable_typo'] += 1
                        measurements_df.at[index, CalculatedColumns.FLAG.value] += SewageFlag.SEWAGE_FLOW_PROBABLE_TYPO.value
                    # is the mean flow a factor of N (default: 2) higher than the dry weather flow --> heavy precipitation
                    elif (current_mean_flow / dry_flow) > self.heavy_precipitation_factor:
                        mean_flow_quality_stats['heavy_precipitation'] += 1
                        measurements_df.at[index, CalculatedColumns.FLAG.value] += SewageFlag.SEWAGE_FLOW_HEAVY_PRECIPITATION.value
                    # is the mean flow a factor of N (default: 1.5) less than the dry weather flow --> probably typo
                    elif (dry_flow / current_mean_flow) > self.mean_sewage_flow_below_typo_factor:
                        mean_flow_quality_stats['probable_typo'] += 1
                        measurements_df.at[index, CalculatedColumns.FLAG.value] += SewageFlag.SEWAGE_FLOW_PROBABLE_TYPO.value
                    else:
                        mean_flow_quality_stats['passed'] += 1
                else:
                    measurements_df.at[index, CalculatedColumns.FLAG.value] += SewageFlag.SEWAGE_FLOW_NOT_ENOUGH_PREVIOUS_VALUES.value
                    mean_flow_quality_stats['skipped'] += 1
        progress_bar.close()
        self.logger.log.info("[Sewage flow] - [Sample location: '{}'] - "
                             "Mean sewage flow quality control: {}".format(sample_location, mean_flow_quality_stats))
        return dry_flow, is_mean_dry_weather_flow_available

