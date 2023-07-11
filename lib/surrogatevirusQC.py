# Created by alex at 22.06.23
import math
from dateutil.relativedelta import relativedelta
import datetime
from .utils import *
from .plotting import *

class SurrogatevirusQC:

    def __init__(self, interactive, periode_month_surrogatevirus, min_number_surrogatevirus_for_outlier_detection, surrogatevirus_outlier_statistics, output_folder):
        self.interactive = interactive
        self.periode_month_surrogatevirus = periode_month_surrogatevirus
        self.min_number_surrogatevirus_for_outlier_detection = min_number_surrogatevirus_for_outlier_detection
        self.surrogatevirus_outlier_statistics = surrogatevirus_outlier_statistics
        self.output_folder = output_folder
        self.logger = SewageLogger(output_folder)

    def __get_start_timeframe(self, current_date: datetime):
        start_timeframe = (current_date - relativedelta(months=self.periode_month_surrogatevirus)).strftime('%Y-%m-%d')
        return start_timeframe

    def filter_dry_days_time_frame(self, sample_location: str, measurements_df: pd.DataFrame):
        """
         Get rid of rainy days and measurements with empty surrogatevirus fields
         """
        for sVirus in Columns.get_surrogatevirus_columns():
            measurements_df[CalculatedColumns.get_surrogate_flag(sVirus)] += np.where((
                (measurements_df[Columns.TROCKENTAG.value] != "Ja") |
                measurements_df[sVirus].isnull()),
                SewageFlag.SURROGATEVIRUS_VALUE_NOT_USABLE.value, 0)

            self.logger.log.info("[Surrogatevirus: {}] - [Sample location: '{}'] - \n"
                                 "\tNumber of excluded measurements because of rain or no value: {}"
                                 "\tout of {} samples.".
                                 format(sVirus, sample_location, len(measurements_df[measurements_df[CalculatedColumns.get_surrogate_flag(sVirus)] != 0]),
                                        measurements_df.shape[0]))

    def __get_previous_surrogatevirus_values (self, measurements_df: pd.DataFrame, current_measurement, sVirus):
        """
          Timeframe of n month, all surrogatevirus measurements that were not flagged
        """
        #get current timeframe eg. last 4 month from current measurement
        start_timeframe = self.__get_start_timeframe(current_measurement[Columns.DATE.value])

        current_timeframe = measurements_df[(measurements_df[Columns.DATE.value] > start_timeframe) &
                                            (measurements_df[Columns.DATE.value] < current_measurement[Columns.DATE.value])]

        # remove previously flagged values (rain and missing value) and flagged outliers
        surrogate_flag = current_timeframe[CalculatedColumns.get_surrogate_flag(sVirus)]
        surrogate_outlier_flag = current_timeframe[CalculatedColumns.get_surrogate_outlier_flag(sVirus)]
        current_timeframe = current_timeframe[(SewageFlag.is_not_flag_set_for_series(surrogate_flag,SewageFlag.SURROGATEVIRUS_VALUE_NOT_USABLE)) &
                                               (SewageFlag.is_not_flag_set_for_series(surrogate_outlier_flag,SewageFlag.SURROGATEVIRUS_OUTLIER))]
        sVirus_values_to_take = current_timeframe[[Columns.DATE.value, sVirus]]

        return sVirus_values_to_take

    def is_surrogatevirus_outlier(self, sample_location: str, measurements_df: pd.DataFrame):
        """
            Detect surrogatevirus outlier, for each measurement and surogatevirus
        """
        self.logger.log.info(
            "[Surrogatevirus outlier detection] - [Sample location: '{}'] - Using following outlier detection methods: {}".format(
                sample_location, self.surrogatevirus_outlier_statistics))

        surrogatevirus_quality_stats = dict()
        out_list = ['surrogatevirus_total']
        for sVirus in Columns.get_surrogatevirus_columns():
            out_list += [sVirus + '_skipped', sVirus + '_passed', sVirus + '_failed']
        for k in out_list:
            surrogatevirus_quality_stats.setdefault(k, 0)

        progress_bar = self.logger.get_progress_bar(CalculatedColumns.get_num_of_unprocessed(measurements_df), "Analyzing surrogatevirus outliers")
        for index, current_measurement in measurements_df.iterrows():
            if CalculatedColumns.needs_processing(current_measurement):
                progress_bar.update(1)
                surrogatevirus_quality_stats['surrogatevirus_total'] += 1
                for sVirus in Columns.get_surrogatevirus_columns():
                    if SewageFlag.is_not_flag(current_measurement[CalculatedColumns.get_surrogate_flag(sVirus)], SewageFlag.SURROGATEVIRUS_VALUE_NOT_USABLE):
                        sVirus_values_to_take = self.__get_previous_surrogatevirus_values(measurements_df, current_measurement, sVirus)
                        if len(sVirus_values_to_take) > self.min_number_surrogatevirus_for_outlier_detection:
                            is_outlier = detect_outliers(self.surrogatevirus_outlier_statistics, sVirus_values_to_take[sVirus],
                                                 current_measurement[sVirus])
                            if is_outlier:
                                measurements_df.at[index, CalculatedColumns.get_surrogate_outlier_flag(sVirus)] += SewageFlag.SURROGATEVIRUS_OUTLIER.value
                                measurements_df.at[index, CalculatedColumns.FLAG.value] += SewageFlag.SURROGATEVIRUS_OUTLIER.value
                                surrogatevirus_quality_stats[sVirus + '_failed'] += 1
                            else:
                                surrogatevirus_quality_stats[sVirus + '_passed'] += 1
                        else:
                            surrogatevirus_quality_stats[sVirus + '_skipped'] += 1
                    else:
                        surrogatevirus_quality_stats[sVirus + '_skipped'] += 1
        progress_bar.close()
        self.logger.log.info(
                "[Surrogatevirus outlier detection] - [Sample location: '{}']\n Final outcome: {}".format(
                    sample_location, surrogatevirus_quality_stats))

        plot_surrogatvirus(measurements_df, sample_location, os.path.join(self.output_folder, "plots", "surrogatvirus"),
                           self.surrogatevirus_outlier_statistics, self.interactive)

        #TODO Datensatz aussortieren, wenn beide Surrogatviren geflagged



