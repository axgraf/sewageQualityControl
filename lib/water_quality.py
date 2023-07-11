# Created by alex at 29.06.23
import pandas as pd
from .constant import *
import math
from .utils import *
from .plotting import *


class WaterQuality:
    def __init__(self, output_folder, water_quality_number_of_last_month,
                 min_number_of_last_measurements_for_water_qc, water_qc_outlier_statistics):
        self.output_folder = output_folder
        self.water_quality_number_of_last_month = water_quality_number_of_last_month
        self.min_number_of_last_measurements_for_water_qc = min_number_of_last_measurements_for_water_qc
        self.water_qc_outlier_statistics = water_qc_outlier_statistics
        self.logger = SewageLogger(self.output_folder)

    def check_water_quality(self, sample_location, measurements_df: pd.DataFrame):
        self.__detect_outliers_in_ammonium(sample_location, measurements_df)
        self.__detect_outliers_in_conductivity(sample_location, measurements_df)
        plot_water_quality(measurements_df, sample_location, os.path.join(self.output_folder, "plots", "water_quality"), self.water_qc_outlier_statistics)

    def __detect_outliers_in_ammonium(self, sample_location, measurements_df: pd.DataFrame):
        self.logger.log.info(
            "[Water quality control: Ammonium] - [Sample location: '{}'] - Using following outlier detection methods: {}".format(
                sample_location, self.water_qc_outlier_statistics))
        water_quality_stats = dict()
        for k in ['skipped', 'passed', 'failed', 'total']:
            water_quality_stats.setdefault(k, 0)
        progress_bar = self.logger.get_progress_bar(CalculatedColumns.get_num_of_unprocessed(measurements_df), "Analyzing ammonium outliers")
        for index, current_measurement in measurements_df.iterrows():
            if CalculatedColumns.needs_processing(current_measurement):
                progress_bar.update(1)
                water_quality_stats['total'] += 1
                last_values = get_last_N_month_and_days(measurements_df, current_measurement, Columns.AMMONIUM.value,
                                               self.water_quality_number_of_last_month, 0, SewageFlag.AMMONIUM_OUTLIER)
                enough_last_values = last_values.shape[0] >= self.min_number_of_last_measurements_for_water_qc
                if current_measurement[Columns.AMMONIUM.value] and not math.isnan(current_measurement[Columns.AMMONIUM.value]):  # only if current value is not empty
                    if not enough_last_values:
                        water_quality_stats['skipped'] += 1
                        measurements_df.at[index, Columns.FLAG.value] += SewageFlag.NOT_ENOUGH_AMMONIUM_VALUES.value
                    else:
                        is_outlier = detect_outliers(self.water_qc_outlier_statistics, last_values[Columns.AMMONIUM.value], current_measurement[Columns.AMMONIUM.value])
                        if is_outlier:
                            measurements_df.at[index, Columns.FLAG.value] += SewageFlag.AMMONIUM_OUTLIER.value
                            water_quality_stats['failed'] += 1
                        else:
                            water_quality_stats['passed'] += 1
                else:
                    water_quality_stats['skipped'] += 1
                    self.logger.log.debug(
                        "[Water quality control: Ammonium] - [Sample location: '{}'] - [Collection date: '{}'] - "
                        "Empty ammonium value: '{}'. Biomarker ratio will be skipped."
                        .format(sample_location, current_measurement[Columns.DATE.value],
                                current_measurement[Columns.AMMONIUM.value]))
        self.logger.log.info("[Water quality control: Ammonium] - [Sample location: '{}'] - "
                             "Ammonium quality control: {}".format(sample_location, water_quality_stats))

    def __detect_outliers_in_conductivity(self, sample_location, measurements_df: pd.DataFrame):
        self.logger.log.info(
            "[Water quality control: Conductivity] - [Sample location: '{}'] - Using following outlier detection methods: {}".format(
                sample_location, self.water_qc_outlier_statistics))
        water_quality_stats = dict()
        for k in ['skipped', 'passed', 'failed', 'total']:
            water_quality_stats.setdefault(k, 0)
        progress_bar = self.logger.get_progress_bar(CalculatedColumns.get_num_of_unprocessed(measurements_df), "Analyzing conductivity outliers")
        for index, current_measurement in measurements_df.iterrows():
            if CalculatedColumns.needs_processing(current_measurement):
                progress_bar.update(1)
                water_quality_stats['total'] += 1
                last_values = get_last_N_month_and_days(measurements_df, current_measurement, Columns.CONDUCTIVITY.value,
                                               self.water_quality_number_of_last_month, 0, SewageFlag.CONDUCTIVITY_OUTLIER)
                enough_last_values = last_values.shape[0] >= self.min_number_of_last_measurements_for_water_qc
                if current_measurement[Columns.CONDUCTIVITY.value] and not math.isnan(current_measurement[Columns.CONDUCTIVITY.value]):  # only if current value is not empty
                    if not enough_last_values:
                        water_quality_stats['skipped'] += 1
                        measurements_df.at[index, CalculatedColumns.FLAG.value] += SewageFlag.NOT_ENOUGH_CONDUCTIVITY_VALUES.value
                    else:
                        is_outlier = detect_outliers(self.water_qc_outlier_statistics, last_values[Columns.CONDUCTIVITY.value], current_measurement[Columns.CONDUCTIVITY.value])
                        if is_outlier:
                            measurements_df.at[index, CalculatedColumns.FLAG.value] += SewageFlag.CONDUCTIVITY_OUTLIER.value
                            water_quality_stats['failed'] += 1
                        else:
                            water_quality_stats['passed'] += 1
                else:
                    water_quality_stats['skipped'] += 1
                    self.logger.log.debug(
                        "[Water quality control: Conductivity] - [Sample location: '{}'] - [Collection date: '{}'] - "
                        "Empty conductivity value: '{}'. Biomarker ratio will be skipped."
                        .format(sample_location, current_measurement[Columns.DATE.value],
                                current_measurement[Columns.CONDUCTIVITY.value]))
        progress_bar.close()
        self.logger.log.info(
            "[Water quality control: Conductivity] - [Sample location: '{}'] - Conductivity quality control: {}".format(
                sample_location, water_quality_stats))

