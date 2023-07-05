# Created by alex at 15.06.23

from enum import Enum, Flag, auto

import pandas as pd


class SewageFlag(Flag):
    ### general flags ###
    COMMENT_NOT_EMPTY = auto()

    ### biomarker flags ###
    MIN_BIOMARKER_NUMBER_NOT_REACHED = auto()
    BIOMARKER_BELOW_THRESHOLD = auto()
    BIOMARKER_OUTLIER = auto()

    ### biomarker ratio flags ###
    NOT_ENOUGH_PREVIOUS_BIOMARKER_VALUES = auto()
    BIOMARKER_RATIO_OUTLIER = auto()

    ### surrogatvirus flags ###
    SURROGATEVIRUS_VALUE_NOT_USABLE = auto()
    SURROGATEVIRUS_OUTLIER = auto()

    ### water quality ####
    NOT_ENOUGH_AMMONIUM_VALUES = auto()
    AMMONIUM_OUTLIER = auto()
    NOT_ENOUGH_CONDUCTIVITY_VALUES = auto()
    CONDUCTIVITY_OUTLIER = auto()

    @staticmethod
    def is_not_flag_set_for_series(series: pd.Series, sewage_flag) -> pd.Series:
        """
            checks if the given flag is NOT contained in a series/dataframe and returns a boolean series
        """
        if isinstance(sewage_flag, SewageFlag):
            return (series & sewage_flag.value) == 0
        else:
            raise ValueError("The given flag is not of type 'SewageFlag'")

    @staticmethod
    def is_flag_set_for_series(series: pd.Series, sewage_flag) -> pd.Series:
        """
            checks if the given flag is contained in the series and returns a boolean series
        """
        if isinstance(sewage_flag, SewageFlag):
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
    CROSSPHAGE = "crassphage"
    PMMOV = "pmmov"
    TROCKENTAG = "trockentag"

    # additional columns for analyses #
    FLAG = "flag"
    BIOMARKER_FLAG = "flag_biomaker"
    BIOMARKER_RATIO_FLAG = "flag_ratio"
    NUMBER_OF_USABLE_BIOMARKERS = "num_usable_biomarkers"
    SURROGATEVIRUS_FLAG = "flag_surrogatevirus"
    SURROGATEVIRUS_OUTLIER_FLAG = "flag_surrogatevirus_outlier"

    @staticmethod
    def get_biomarker_columns():
        biomarker_columns = [Columns.BIOMARKER_N1.value, Columns.BIOMARKER_N2.value,
                             Columns.BIOMARKER_N3.value, Columns.BIOMARKER_E.value,
                             Columns.BIOMARKER_ORF.value, Columns.BIOMARKER_RDRP.value]
        return biomarker_columns

    @staticmethod
    def get_biomarker_flag(biomarker):
        return Columns.BIOMARKER_FLAG.value + "_" + biomarker

    @staticmethod
    def get_biomaker_ratio_flag(biomarker1, biomarker2):
        return Columns.BIOMARKER_RATIO_FLAG.value + "_" + biomarker1 + "/" + biomarker2

    @staticmethod
    def get_surrogatevirus_columns():
        surrogatevirus_columns = [Columns.CROSSPHAGE.value, Columns.PMMOV.value]
        return surrogatevirus_columns
    @staticmethod
    def get_surrogate_outlier_flag(sVirus):
        return Columns.SURROGATEVIRUS_OUTLIER_FLAG.value + "_" + sVirus

    @staticmethod
    def get_surrogate_flag(sVirus):
        return Columns.SURROGATEVIRUS_FLAG.value + "_" + sVirus