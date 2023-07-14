# Created by alex at 29.06.23
import pandas as pd
from .constant import *
import math
from .utils import *
from .statistics import *
from .plotting import *


class WaterQuality:
    def __init__(self, output_folder, sewageStat: SewageStat, water_quality_number_of_last_month,
                 min_number_of_last_measurements_for_water_qc, water_qc_outlier_statistics):
        self.output_folder = output_folder
        self.sewageStat = sewageStat
        self.water_quality_number_of_last_month = water_quality_number_of_last_month
        self.min_number_of_last_measurements_for_water_qc = min_number_of_last_measurements_for_water_qc
        self.water_qc_outlier_statistics = water_qc_outlier_statistics
        self.logger = SewageLogger(self.output_folder)

    def check_water_quality(self, sample_location, measurements: pd.DataFrame, index):
        self.__detect_outliers_in_ammonium(sample_location, measurements, index)
        self.__detect_outliers_in_conductivity(sample_location, measurements, index)

    def __detect_outliers_in_ammonium(self, sample_location, measurements: pd.DataFrame, index):
        current_measurement = measurements.iloc[index]
        last_values = get_last_N_month_and_days(measurements, current_measurement, Columns.AMMONIUM.value,
                                       self.water_quality_number_of_last_month, 0, SewageFlag.AMMONIUM_OUTLIER)
        enough_last_values = last_values.shape[0] >= self.min_number_of_last_measurements_for_water_qc
        if current_measurement[Columns.AMMONIUM.value] and not math.isnan(current_measurement[Columns.AMMONIUM.value]):  # only if current value is not empty
            if not enough_last_values:
                self.sewageStat.add_water_quality_outlier("Ammonium", 'skipped')
                measurements.at[index, CalculatedColumns.FLAG.value] += SewageFlag.NOT_ENOUGH_AMMONIUM_VALUES.value
            else:
                is_outlier = detect_outliers(self.water_qc_outlier_statistics, last_values[Columns.AMMONIUM.value], current_measurement[Columns.AMMONIUM.value])
                if is_outlier:
                    measurements.at[index, CalculatedColumns.FLAG.value] += SewageFlag.AMMONIUM_OUTLIER.value
                    self.sewageStat.add_water_quality_outlier("Ammonium", 'failed')
                else:
                    self.sewageStat.add_water_quality_outlier("Ammonium", 'passed')
        else:
            self.sewageStat.add_water_quality_outlier("Ammonium", 'skipped')


    def __detect_outliers_in_conductivity(self, sample_location, measurements: pd.DataFrame, index):
        current_measurement = measurements.iloc[index]
        last_values = get_last_N_month_and_days(measurements, current_measurement, Columns.CONDUCTIVITY.value,
                                       self.water_quality_number_of_last_month, 0, SewageFlag.CONDUCTIVITY_OUTLIER)
        enough_last_values = last_values.shape[0] >= self.min_number_of_last_measurements_for_water_qc
        if current_measurement[Columns.CONDUCTIVITY.value] and not math.isnan(current_measurement[Columns.CONDUCTIVITY.value]):  # only if current value is not empty
            if not enough_last_values:
                self.sewageStat.add_water_quality_outlier("Conductivity", 'skipped')
                measurements.at[index, CalculatedColumns.FLAG.value] += SewageFlag.NOT_ENOUGH_CONDUCTIVITY_VALUES.value
            else:
                is_outlier = detect_outliers(self.water_qc_outlier_statistics, last_values[Columns.CONDUCTIVITY.value], current_measurement[Columns.CONDUCTIVITY.value])
                if is_outlier:
                    measurements.at[index, CalculatedColumns.FLAG.value] += SewageFlag.CONDUCTIVITY_OUTLIER.value
                    self.sewageStat.add_water_quality_outlier("Conductivity", 'failed')
                else:
                    self.sewageStat.add_water_quality_outlier("Conductivity", 'passed')
        else:
            self.sewageStat.add_water_quality_outlier("Conductivity", 'skipped')

