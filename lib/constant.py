# Created by alex at 15.06.23

from enum import Enum, Flag

import numpy as np
import pandas as pd


class SewageFlag(Flag):
    ### general flags ###
    COMMENT_NOT_EMPTY = 1
    MISSING_MEAN_SEWAGE_FLOW = 2

    ### biomarker flags ###
    BIOMARKER_BELOW_THRESHOLD_OR_EMPTY = 4
    BIOMARKER_PROBABLE_OUTLIER = 8
    BIOMARKER_VALIDATED_OUTLIER = 16

    ### biomarker ratio flags ###
    NOT_ENOUGH_PREVIOUS_BIOMARKER_VALUES = 32
    BIOMARKER_RATIO_OUTLIER = 64

    ### flow rate flags ###
    SEWAGE_FLOW_NOT_ENOUGH_PREVIOUS_VALUES = 128
    SEWAGE_FLOW_HEAVY_PRECIPITATION = 256
    SEWAGE_FLOW_PROBABLE_TYPO = 512

    ### surrogatvirus flags ###
    SURROGATEVIRUS_VALUE_NOT_USABLE = 1024
    SURROGATEVIRUS_OUTLIER = 2048

    ### water quality ####
    NOT_ENOUGH_AMMONIUM_VALUES = 4096
    AMMONIUM_OUTLIER = 8192
    NOT_ENOUGH_CONDUCTIVITY_VALUES = 16384
    CONDUCTIVITY_OUTLIER = 32768

    ### biomarker normalization ###
    NOT_ENOUGH_BIOMARKERS_FOR_NORMALIZATION = 65536
    REPRODUCTION_NUMBER_OUTLIER_SKIPPED = 131072
    REPRODUCTION_NUMBER_OUTLIER = 262144

    def __str__(self):
        return self.name

    def __eq__(self, other):
        return self.name == other.name

    @staticmethod
    def remove_flag_from_index_column(df, index, flag_column, sewage_flag) -> None:
        if isinstance(sewage_flag, SewageFlag):
            current_flag = df.iloc[index][flag_column]
            has_flag = SewageFlag.is_flag(current_flag, sewage_flag)
            if has_flag:
                df.at[index, flag_column] -= sewage_flag.value
        else:
            raise ValueError("The given flag is not of type 'SewageFlag'")

    @staticmethod
    def add_flag_to_index_column(df, index, flag_column, sewage_flag) -> None:
        if isinstance(sewage_flag, SewageFlag):
            current_flag = df.iloc[index][flag_column]
            has_flag = SewageFlag.is_flag(current_flag, sewage_flag)
            if not has_flag:
                df.at[index, flag_column] += sewage_flag.value
        else:
            raise ValueError("The given flag is not of type 'SewageFlag'")

    @staticmethod
    def add_series_flag_to_column(df: pd.DataFrame, flag_column, sewage_flag_series: np.array) -> None:
        for index, row in df.iterrows():
            current_flag = row[flag_column]
            new_flag = SewageFlag.get_flag_from_value(sewage_flag_series[index])
            if not SewageFlag.is_flag(current_flag, new_flag):
                df.at[index, flag_column] += new_flag.value

    @staticmethod
    def is_not_flag_set_for_series(series: pd.Series, sewage_flag) -> pd.Series:
        """
            checks if the given flag is NOT contained in a series/dataframe and returns a boolean series
        """
        if isinstance(sewage_flag, SewageFlag):
            series = series.astype(int)
            return (series & sewage_flag.value) == 0
        else:
            raise ValueError("The given flag is not of type 'SewageFlag'")

    @staticmethod
    def is_flag_set_for_series(series: pd.Series, sewage_flag) -> pd.Series:
        """
            checks if the given flag is contained in the series and returns a boolean series
        """
        if isinstance(sewage_flag, SewageFlag):
            series = series.astype(int)
            return (series & sewage_flag.value) == sewage_flag.value
        else:
            raise ValueError("The given flag is not of type 'SewageFlag'")

    @staticmethod
    def explain_flag_series(series: pd.Series) -> pd.Series:
        """
        Convert a flag column to named flags.
        E.g.: SewageFlag.convert_flag_series_to_named(measurements[Columns.FLAG.value])
        to convert the integer flag column to named flag column
        """
        return series.apply(SewageFlag.explain_flags_to_string)

    @staticmethod
    def explain_flags_to_array(flag_value: int) -> []:
        """
        Returns a list with names of all flags that are stored in the given integer flag.
        """
        named_flags = []
        for f in SewageFlag:
            if f.value & flag_value == f.value:
                named_flags.append(f.name)
        return named_flags

    @staticmethod
    def explain_flags_to_string(flag_value: int) -> str:
        """
        Returns a string representation of all flags that are stored in the given flag_value.
        """
        named_flags = SewageFlag.explain_flags_to_array(flag_value)
        return ', '.join(named_flags)

    @staticmethod
    def is_flag(flag_value: int, sewage_flag) -> bool:
        """
        Check if a flag is set
        """
        if isinstance(sewage_flag, SewageFlag):
            return (sewage_flag.value & flag_value) == sewage_flag.value
        else:
            raise ValueError("The given flag is not of type 'SewageFlag'")

    @staticmethod
    def is_not_flag(flag_value: int, sewage_flag) -> bool:
        """
        Check if a flag is NOT set
        """
        if isinstance(sewage_flag, SewageFlag):
            return (sewage_flag.value & flag_value) == 0
        else:
            raise ValueError("The given flag is not of type 'SewageFlag'")

    @staticmethod
    def get_flag_from_value(flag_value: int):
        return SewageFlag(flag_value)


class CalculatedColumns(Enum):
    # calculated columns for analyses #
    FLAG = "flag", np.int
    BIOMARKER_FLAG = "flag_biomaker", np.int
    BIOMARKER_RATIO_FLAG = "flag_ratio", np.int
    NUMBER_OF_USABLE_BIOMARKERS = "num_usable_biomarkers", np.int
    SURROGATEVIRUS_FLAG = "flag_surrogatevirus", np.int
    SURROGATEVIRUS_OUTLIER_FLAG = "flag_surrogatevirus_outlier", np.int
    NORMALIZED_MEAN_BIOMARKERS = "normalized_mean_biomarkers", np.float
    BASE_REPRODUCTION_FACTOR = "base_reproduction_factor", np.float
    NEEDS_PROCESSING = "needs_processing", np.bool
    USABLE = "usable", np.bool
    OUTLIER_REASON = "outlier_reason", str

    def __init__(self, _: str, dtype):
        self._type_ = dtype

    def __str__(self):
        return self.value[0]

    @property
    def type(self):
        return self._type_

    @property
    def value(self):
        return self._value_[0]
    @staticmethod
    def get_biomarker_flag_columns() -> []:
        biomarker_flag_columns = [CalculatedColumns.get_biomarker_flag(Columns.BIOMARKER_N1.value),
                                  CalculatedColumns.get_biomarker_flag(Columns.BIOMARKER_N2.value),
                                  CalculatedColumns.get_biomarker_flag(Columns.BIOMARKER_N3.value),
                                  CalculatedColumns.get_biomarker_flag(Columns.BIOMARKER_E.value),
                                  CalculatedColumns.get_biomarker_flag(Columns.BIOMARKER_ORF.value),
                                  CalculatedColumns.get_biomarker_flag(Columns.BIOMARKER_RDRP.value)]
        return biomarker_flag_columns

    @staticmethod
    def needs_processing(row: pd.Series) -> bool:
        if CalculatedColumns.NEEDS_PROCESSING.value in row:
            return row[CalculatedColumns.NEEDS_PROCESSING.value]
        else:
            return True

    @staticmethod
    def get_num_of_unprocessed(measurements_df: pd.DataFrame) -> int:
        if CalculatedColumns.NEEDS_PROCESSING.value in measurements_df:
            needs_processing_df = measurements_df[measurements_df[CalculatedColumns.NEEDS_PROCESSING.value] == True]
            return needs_processing_df.shape[0]
        else:
            return measurements_df.shape[0]
    @staticmethod
    def get_biomarker_flag(biomarker) -> str:
        return CalculatedColumns.BIOMARKER_FLAG.value + "_" + biomarker

    @staticmethod
    def get_biomaker_ratio_flag(biomarker1, biomarker2) -> str:
        return CalculatedColumns.BIOMARKER_RATIO_FLAG.value + "_" + biomarker1 + "/" + biomarker2

    @staticmethod
    def get_surrogate_outlier_flag(sVirus) -> str:
        return CalculatedColumns.SURROGATEVIRUS_OUTLIER_FLAG.value + "_" + sVirus

    @staticmethod
    def get_surrogate_flag(sVirus) -> str:
        return CalculatedColumns.SURROGATEVIRUS_FLAG.value + "_" + sVirus

    @staticmethod
    def add_outlier_reason(measurements_df: pd.DataFrame, index: int, reason: str) -> None:
        values = measurements_df.iloc[index][CalculatedColumns.OUTLIER_REASON.value] + ", " + reason
        if values.startswith(','):
            values = values[2:]
        measurements_df.at[index, CalculatedColumns.OUTLIER_REASON.value] = values


class Columns(Enum):
    # columns used in sample dataframe #
    DATE = "collectionDate"
    COMMENT_ANALYSIS = "bem_lab"
    COMMENT_OPERATION = "bem_pn"
    BIOMARKER_N1 = 'biomarker_N1'
    BIOMARKER_N2 = 'biomarker_N2'
    BIOMARKER_N3 = 'biomarker_N3'
    BIOMARKER_E = 'biomarker_E'
    BIOMARKER_ORF = 'biomarker_ORF'
    BIOMARKER_RDRP = 'biomarker_RDRP'
    AMMONIUM = "nh4n"
    CONDUCTIVITY = "lf"
    MEAN_SEWAGE_FLOW = "mean_sewage_flow"
    CROSSPHAGE = "crassphage"
    PMMOV = "pmmov"
    TROCKENTAG = "trockentag"

    @staticmethod
    def get_biomarker_columns():
        biomarker_columns = [Columns.BIOMARKER_N1.value, Columns.BIOMARKER_N2.value,
                             Columns.BIOMARKER_N3.value, Columns.BIOMARKER_E.value,
                             Columns.BIOMARKER_ORF.value, Columns.BIOMARKER_RDRP.value]
        return biomarker_columns

    @staticmethod
    def get_surrogatevirus_columns():
        surrogatevirus_columns = [Columns.CROSSPHAGE.value, Columns.PMMOV.value]
        return surrogatevirus_columns
   
