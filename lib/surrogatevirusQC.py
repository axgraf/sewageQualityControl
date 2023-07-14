# Created by alex at 22.06.23
import math
from dateutil.relativedelta import relativedelta
from .utils import *
from .statistics import *
from .plotting import *


class SurrogateVirusQC:

    def __init__(self, sewageStat: SewageStat, periode_month_surrogatevirus, min_number_surrogatevirus_for_outlier_detection, surrogatevirus_outlier_statistics, output_folder):
        self.sewageStat = sewageStat
        self.periode_month_surrogatevirus = periode_month_surrogatevirus
        self.min_number_surrogatevirus_for_outlier_detection = min_number_surrogatevirus_for_outlier_detection
        self.surrogatevirus_outlier_statistics = surrogatevirus_outlier_statistics
        self.output_folder = output_folder
        self.logger = SewageLogger(output_folder)

    def __get_start_timeframe(self, current_date: datetime):
        start_timeframe = (current_date - relativedelta(months=self.periode_month_surrogatevirus)).strftime('%Y-%m-%d')
        return start_timeframe

    def filter_dry_days_time_frame(self, sample_location: str, measurements: pd.DataFrame, index):
        """
         Get rid of rainy days and measurements with empty surrogatevirus fields
         """
        current_measurement = measurements.iloc[index]
        for sVirus in Columns.get_surrogatevirus_columns():
            if current_measurement[Columns.TROCKENTAG.value].lower().strip() != "ja" or math.isnan(current_measurement[sVirus]):
                SewageFlag.add_flag_to_index_column(measurements, index, CalculatedColumns.get_surrogate_flag(sVirus),
                                                    SewageFlag.SURROGATEVIRUS_VALUE_NOT_USABLE)

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

    def is_surrogatevirus_outlier(self, sample_location: str, measurements: pd.DataFrame, index):
        """
            Detect surrogatevirus outlier, for each measurement and surogatevirus
        """
        current_measurement = measurements.iloc[index]
        for sVirus in Columns.get_surrogatevirus_columns():
            if SewageFlag.is_not_flag(current_measurement[CalculatedColumns.get_surrogate_flag(sVirus)], SewageFlag.SURROGATEVIRUS_VALUE_NOT_USABLE):
                sVirus_values_to_take = self.__get_previous_surrogatevirus_values(measurements, current_measurement, sVirus)
                if len(sVirus_values_to_take) > self.min_number_surrogatevirus_for_outlier_detection:
                    is_outlier = detect_outliers(self.surrogatevirus_outlier_statistics, sVirus_values_to_take[sVirus],
                                                 current_measurement[sVirus])
                    if is_outlier:
                        measurements.at[index, CalculatedColumns.get_surrogate_outlier_flag(sVirus)] += SewageFlag.SURROGATEVIRUS_OUTLIER.value
                        measurements.at[index, CalculatedColumns.FLAG.value] += SewageFlag.SURROGATEVIRUS_OUTLIER.value
                        self.sewageStat.add_surrogate_virus_outlier(sVirus, 'outlier')
                    else:
                        self.sewageStat.add_surrogate_virus_outlier(sVirus, 'passed')

                else:
                    self.sewageStat.add_surrogate_virus_outlier(sVirus, 'skipped')

        #TODO Datensatz aussortieren, wenn beide Surrogatviren geflagged



