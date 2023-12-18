# Created by alex at 20.06.23
import os.path
import sys
from typing import List
import datetime
import dateutil.relativedelta
import numpy as np
import pandas as pd
import logging
import tqdm
import scipy.stats as st
from sklearn.ensemble import IsolationForest
from sklearn.neighbors import LocalOutlierFactor
from sklearn.svm import OneClassSVM
from .sewage import SewageSample
from .constant import *

column_map = {
    'ANFANG': 'collectionDate',
    'BEM_LAB': 'bem_lab',
    'BEM_PN': 'bem_pn',
    'N1_LAB': 'biomarker_N1',
    'N2_LAB': 'biomarker_N2',
    'N3_LAB': 'biomarker_N3',
    'E_LAB': 'biomarker_E',
    'ORF_LAB': 'biomarker_ORF',
    'RDRP_LAB': 'biomarker_RDRP',
    'NH4N': 'nh4n',
    'LF': "lf",
    'VOLUMENSTROM': 'mean_sewage_flow',
    'CRASSPHAGE': 'crassphage',
    'PMMOV': 'pmmov',
    'TRO_TAG': 'trockentag'
}


def read_excel_input_files(input_file: str):
    f = pd.ExcelFile(input_file)
    location_dict = dict()
    for sheet in f.sheet_names:
        df = f.parse(sheet)
        df.rename(columns=column_map, inplace=True)
        df = df[list(column_map.values())]
        location_dict[sheet] = df
    return location_dict


def convert_sample_list2pandas(measurements: List[SewageSample]):
    table = pd.DataFrame.from_records([vars(s) for s in measurements])
    return table


def detect_outliers(outlier_statistics, train_values, test_value, isFactor=False):
    outlier_detected = []
    if 'svm' in outlier_statistics or 'all' in outlier_statistics:
        is_svm_outlier = is_oneClassSVM(test_value, train_values)
        outlier_detected.append(is_svm_outlier)
    if 'lof' in outlier_statistics or 'all' in outlier_statistics:
        is_lof_outlier = is_outlier_local_outlier_factor(test_value, train_values)
        outlier_detected.append(is_lof_outlier)
    if 'rf' in outlier_statistics or 'all' in outlier_statistics:
        is_isolation_forest_outlier, isolation_forest_score_ratio = is_outlier_isolation_forest(test_value, train_values)
        outlier_detected.append(is_isolation_forest_outlier)
    if 'ci' in outlier_statistics or 'all' in outlier_statistics:
        is_confidence_interval_99_outlier, confidence_interval = is_confidence_interval_outlier(test_value, train_values, 0.99)
        outlier_detected.append(is_confidence_interval_99_outlier)
    if 'iqr' in outlier_statistics or 'all' in outlier_statistics:
        is_iqr_outlier, iqr_range = interquartile_range(test_value, train_values, isFactor)
        outlier_detected.append(is_iqr_outlier)
    if 'zscore' in outlier_statistics or 'all' in outlier_statistics:
        is_zscore_outlier, z_score_threshold = is_outlier_modified_z_score(test_value, train_values)
        outlier_detected.append(is_zscore_outlier)
    return all(outlier_detected)



def is_outlier_local_outlier_factor(test_value: float, train_values: List[float], contamination='auto'):
    """
    The Local Outlier Factor measures the local deviation of the density of a given sample with respect to its neighbors.
    It is local in that the anomaly score depends on how isolated the object is with respect to the surrounding neighborhood.
    More precisely, locality is given by k-nearest neighbors, whose distance is used to estimate the local density.
    By comparing the local density of a sample to the local densities of its neighbors, one can identify samples
    that have a substantially lower density than their neighbors. These are considered outliers.

    :param test_value: the value to be checked as an outlier
    :param train_values: values to be used for training the model
    :param contamination:  the amount of contamination, i.e. the proportion of outliers in the data set. Range should be: (0, 0.5]
    """
    X = np.array(train_values).reshape(-1, 1)
    neighbours = 20 if len(X) > 20 else len(X) - 1
    lof_novelty = LocalOutlierFactor(n_neighbors=neighbours, novelty=True, contamination=contamination,).fit(X)
    test_value = np.array(test_value).reshape(1, -1)
    prediction = lof_novelty.predict(test_value)
    return prediction[0] == -1


def is_oneClassSVM(test_value: float, train_values: List[float]):
    """
    One-class SVM with non-linear kernel (RBF). One-class SVM is an unsupervised algorithm
    that learns a decision function for novelty detection: classifying new data as similar or different to the training set.
    :param test_value: the value to be checked as an outlier
    :param train_values: values to be used for training the model
    """
    model = OneClassSVM(nu=0.1, kernel="rbf", gamma=0.2)
    model.fit(np.array(train_values).reshape(-1, 1))
    prediction = model.predict(np.array(test_value).reshape(-1, 1))   #1 = inlier;  -1 : outlier
    print(prediction)
    if prediction[0] == 1:
        return False
    return True

def is_outlier_isolation_forest(test_value: float, train_values: List[float], contamination=0.1):
    """
    The IsolationForest ‘isolates’ observations by randomly selecting a feature and then randomly selecting a
    split value between the maximum and minimum values of the selected feature. Since recursive partitioning can be
    represented by a tree structure, the number of splittings required to isolate a sample is equivalent to the
    path length from the root node to the terminating node. This path length, averaged over a forest
    of such random trees, is a measure of normality and our decision function. Random partitioning produces noticeably
    shorter paths for anomalies. Hence, when a forest of random trees collectively produce shorter path lengths
    for particular samples, they are highly likely to be anomalies.

    :param test_value: the value to be checked as an outlier
    :param train_values: values to be used for training the model
    :param contamination:  the amount of contamination, i.e. the proportion of outliers in the data set. Range should be: (0, 0.5]
    """
    X = pd.DataFrame(train_values)
    X.rename(columns={X.columns[0]: 'samples'}, inplace=True)
    model = IsolationForest(n_estimators=100, warm_start=False, contamination=contamination,
                            n_jobs=4)  # contamination="auto" else range should be (0, 0.5]
    model.fit(X.values)
    test = np.array(test_value).reshape(1, -1)
    score = model.decision_function(test)
    outlier = model.predict(test)
    return outlier[0] == -1, score


def is_outlier_modified_z_score(test_value: float, train_values: List[float]):
    train_values = train_values.tolist()
    median = np.median(train_values)
    abs_diff = np.abs(train_values - median)
    median_abs_diff = np.median(abs_diff)
    modified_zscore = 0.6745 * ((test_value - median) / median_abs_diff)
    max_standard_deviation = 3.5  # how many std deviations; good value
    return modified_zscore > max_standard_deviation, modified_zscore


def interquartile_range(test_value: float, train_values: List[float], isFactor=False):
    train_values = train_values.tolist()
    q1 = np.quantile(train_values, 0.25)
    q3 = np.quantile(train_values, 0.75)
    iqr = q3 - q1
    # calculation for multiplier selection based on the standard deviation
    # std_dev = 3.5
    # multiplier = (std_dev - 0.675)/(0.675 + 0.675)
    # e.g. IQR multiplier of 1.7 = standard deviation of 3; 1.5 = standard deviation of 2.7
    multiplier = 1.5
    minimum = q1 - multiplier * iqr
    if isFactor:
        minimum = 1 / (q3 + multiplier * iqr)
    maximum = q3 + multiplier * iqr
    if minimum <= test_value <= maximum:
        return False, (minimum, maximum)
    return True, (minimum, maximum)


def is_confidence_interval_outlier(test_value: float, train_values: List[float], confidence: float):
    train_values = train_values.tolist()
    if len(train_values) < 30:  # use t-distribution in case less than 30 samples are used
        confidence_interval = st.t.interval(alpha=confidence, df=len(train_values) - 1,
                                            loc=np.mean(train_values),
                                            scale=st.sem(train_values))
    else:  # use normal distribution
        confidence_interval = st.norm.interval(alpha=confidence,
                                               loc=np.mean(train_values),
                                               scale=st.sem(train_values))
    if confidence_interval[0] <= test_value <= confidence_interval[1]:
        return False, confidence_interval
    return True, confidence_interval


def get_last_values(measurements_df: pd.DataFrame, index, column_name,
                              sewage_flag: SewageFlag = None, additional_sewage_flag: SewageFlag = None):
    """
    Obtain last values from last N month from the data frame. In case a flag is provided values that do have the flag set
    will be filtered.

    :param measurements_df: full data frame ot select last values
    :param index of current measurement
    :param column_name: remove NA from selected column
    :param sewage_flag: Sewage flag to filter for previous outliers; must be not set in flags
    :param additional_sewage_flag: Additional Sewage flag to filter for previous outliers; must be not set in flags
    :return: data frame with entries from last N month
    """
    max_idx = (index - 1 - self.max_number_biomarkers_for_outlier_detection) \
        if (index - 1 - self.max_number_biomarkers_for_outlier_detection) >= 0 else 0
    min_idx = index if index >= 0 else 0
    last_values = measurements_df.iloc[max_idx: min_idx]
    # remove outliers based on sewage_flag -> filters values which do not have the flag
    if sewage_flag:
        last_values = last_values[SewageFlag.is_not_flag_set_for_series(last_values[CalculatedColumns.FLAG.value], sewage_flag)]
    if additional_sewage_flag:
        last_values = last_values[SewageFlag.is_not_flag_set_for_series(last_values[CalculatedColumns.FLAG.value], additional_sewage_flag)]
    last_values = last_values[last_values[column_name].notna()]
    return last_values


def get_last_N_month_and_days(measurements_df: pd.DataFrame, current_measurement, column_name, num_month, num_days,
                              sewage_flag: SewageFlag = None, additional_sewage_flag: SewageFlag = None):
    """
    Obtain last values from last N month from the data frame. In case a flag is provided values that do have the flag set
    will be filtered.

    :param measurements_df: full data frame ot select last values
    :param current_measurement: current value
    :param column_name: remove NA from selected column
    :param num_month: number of last month to select values
    :param num_days: number of last days to select values
    :param sewage_flag: Sewage flag to filter for previous outliers; must be not set in flags
    :param additional_sewage_flag: Additional Sewage flag to filter for previous outliers; must be not set in flags
    :return: data frame with entries from last N month
    """
    min_date = current_measurement[Columns.DATE.value]
    if num_month > 0:
        min_date = min_date + dateutil.relativedelta.relativedelta(months=-num_month)
    if num_days > 0:
        min_date = min_date + dateutil.relativedelta.relativedelta(days=-num_days)
    last_values = measurements_df[(measurements_df[Columns.DATE.value] >= min_date) &
                                  (measurements_df[Columns.DATE.value] < current_measurement[Columns.DATE.value])]
    # remove outliers based on sewage_flag -> filters values which do not have the flag
    if sewage_flag:
        last_values = last_values[SewageFlag.is_not_flag_set_for_series(last_values[CalculatedColumns.FLAG.value], sewage_flag)]
    if additional_sewage_flag:
        last_values = last_values[SewageFlag.is_not_flag_set_for_series(last_values[CalculatedColumns.FLAG.value], additional_sewage_flag)]
    last_values = last_values[last_values[column_name].notna()]
    return last_values


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


class SewageLogger:
    _instance = None
    log = None
    verbosity = None
    quiet = None

    def __init__(self, output_folder, verbosity=None, quiet=False):
        self.output_folder = os.path.join(output_folder, "logs")
        self.verbosity = verbosity
        self.quiet = quiet
        self.__initalize()

    def __initalize(self):
        if not os.path.exists(self.output_folder):
            os.makedirs(self.output_folder)
        if not self.log:
            self.log = self.__setup_logger()

    def __new__(cls, *args, **kwargs):
        if not cls._instance:
            cls._instance = object.__new__(cls)
        return cls._instance

    def logInfo(self, message):

            self.log.info(message)

    def __setup_logger(self):
        logger = logging.getLogger(__name__)
        logger.setLevel(logging.DEBUG)
        fmt = "%(asctime)s - %(levelname)s - %(message)s"
        stdout_handler = logging.StreamHandler()
        stdout_handler.setLevel(logging.INFO)
        if self.verbosity and self.verbosity >= 1:
            stdout_handler.setLevel(logging.DEBUG)
        if self.quiet:
            stdout_handler.setLevel(logging.ERROR)
        stdout_handler.setFormatter(CustomFormatter(fmt))
        today = datetime.date.today()
        log_file = os.path.join(self.output_folder, 'sewage_quality_{}.log'.format(today.strftime('%Y_%m_%d')))
        file_handler = logging.FileHandler(log_file, mode='w')
        file_handler.setLevel(logging.INFO)
        file_handler.setFormatter(logging.Formatter(fmt))
        logger.addHandler(stdout_handler)
        logger.addHandler(file_handler)
        return logger

    def get_progress_bar(self, total, text):
        return tqdm.tqdm(total=total, unit=' samples', colour="blue", ncols=100, desc=text, file=sys.stdout)
