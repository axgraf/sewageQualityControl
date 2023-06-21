# Created by alex at 20.06.23
from .sewage import SewageSample
from typing import List
import numpy as np
import pandas as pd
import logging
import scipy.stats as st


def convert_sample_list2pandas(measurements: List[SewageSample]):
    table = pd.DataFrame.from_records([vars(s) for s in measurements])
    return table


def is_outlier_modified_z_score(test_value: float, values: List[float]):
    median = np.median(values)
    abs_diff = np.abs(values - median)
    median_abs_diff = np.median(abs_diff)
    modified_zscore = 0.6745 * (test_value - median) / median_abs_diff
    return modified_zscore > 3.5, modified_zscore


def inter_quantil_range(test_value: float, values: List[float]):
    q1 = np.quantile(values, 0.25)
    q3 = np.quantile(values, 0.75)
    iqr = q3 - q1
    # IQR multiplier of 1.7 is similar to stdev multiplier of 3, else use 1.5
    minimum = q1 - 1.7 * iqr
    maximum = q3 + 1.7 * iqr
    if minimum <= test_value <= maximum:
        return False, (minimum, maximum)
    return True, (minimum, maximum)


def is_confidence_interval_outlier(test_value: float, values: List[float], confidence: float):
    if len(values) < 30:  ## use t-distribution in case less than 30 samples are used
        confidence_interval = st.t.interval(alpha=confidence, df=len(values) - 1,
                                            loc=np.mean(values),
                                            scale=st.sem(values))
    else:  # use normal distribution
        confidence_interval = st.norm.interval(alpha=confidence,
                                               loc=np.mean(values),
                                               scale=st.sem(values))
    if confidence_interval[0] <= test_value <= confidence_interval[1]:
        return False, confidence_interval
    return True, confidence_interval


def add_default_biomarker_statistic(stat_dict: dict, biomarker1, biomarker2) -> None:
    stat_dict.setdefault(biomarker1 + "/" + biomarker2, dict())
    stat_dict[biomarker1 + "/" + biomarker2].setdefault("passed", 0)
    stat_dict[biomarker1 + "/" + biomarker2].setdefault("failed", 0)
    stat_dict[biomarker1 + "/" + biomarker2].setdefault("skipped", 0)
    stat_dict[biomarker1 + "/" + biomarker2].setdefault("total", 0)


def pretty_print_biomarker_statistic(stat_dict: dict) -> str:
    biomarker_log = ""
    for biomarker, status_dict in stat_dict.items():
        if status_dict['total'] > 0:
            biomarker_log += "\t" + biomarker + ":\t"
            for status, num in status_dict.items():
                biomarker_log += "\t" + status + ":\t" + str(num) + "\t"
            biomarker_log += "\n"
    return biomarker_log

class CustomFormatter(logging.Formatter):
    grey = '\x1b[38;21m'
    green = '\x1b[1;32m'
    blue = '\x1b[38;5;39m'
    yellow = '\x1b[38;5;226m'
    red = '\x1b[38;5;196m'
    bold_red = '\x1b[31;1m'
    purple = '\x1b[1;35m'
    reset = '\x1b[0m'

    def __init__(self, fmt):
        super().__init__()
        self.fmt = fmt
        self.FORMATS = {
            logging.DEBUG: self.grey + self.fmt + self.reset,
            logging.INFO: self.blue + self.fmt + self.reset,
            logging.WARNING: self.yellow + self.fmt + self.reset,
            logging.ERROR: self.red + self.fmt + self.reset,
            logging.CRITICAL: self.bold_red + self.fmt + self.reset
        }

    def format(self, record):
        log_fmt = self.FORMATS.get(record.levelno)
        formatter = logging.Formatter(log_fmt)
        return formatter.format(record)
