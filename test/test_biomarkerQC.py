# Created by alex at 03.07.23
import itertools
import shutil
import datetime
import os.path
from unittest import TestCase
from pandas.testing import *
import numpy as np
import pandas as pd
from lib import biomarkerQC
from lib import constant

test_output_folder = 'tmp'


class TestBiomarkerCommentNotEmpty(TestCase):

    def setUp(self) -> None:
        self.biomarkerQC = biomarkerQC.BiomarkerQC(output_folder=test_output_folder, biomarker_outlier_statistics=None,
                                                   min_number_biomarkers_for_outlier_detection=None,
                                                   min_biomarker_threshold=None, max_number_biomarkers_for_outlier_detection=None, report_number_of_biomarker_outlier=None)
        operation_column = constant.Columns.COMMENT_OPERATION.value
        analysis_column = constant.Columns.COMMENT_ANALYSIS.value
        values = \
            [{operation_column: " ", analysis_column: None},  # no flag
             {operation_column: None, analysis_column: "   "},  # no flag
             {operation_column: None, analysis_column: None},  # no flag
             {operation_column: "wrong method", analysis_column: ""},  # flagged
             {operation_column: "!?# help", analysis_column: None},  # flagged
             {operation_column: None, analysis_column: "analysis went wrong"},  # flagged
             {operation_column: "", analysis_column: "analysis went wrong"},  # flagged
             {operation_column: "222", analysis_column: "analysis went wrong"}]  # flagged

        self.test_df = pd.DataFrame(values)
        self.test_df[constant.Columns.FLAG.value] = 0

    def test_check_comments(self):
        flag_comment_not_empty = constant.SewageFlag.COMMENT_NOT_EMPTY.value
        correct_flags = []
        correct_flags.extend([0] * 3)  # first three should pass
        correct_flags.extend([flag_comment_not_empty] * 5)  # next five should fail
        valid_result = pd.Series(correct_flags, name=constant.Columns.FLAG.value)
        self.biomarkerQC.check_comments('', self.test_df)
        assert_series_equal(self.test_df[constant.Columns.FLAG.value], valid_result)


class TestBiomarkerBelowThreshold(TestCase):
    def setUp(self) -> None:
        self.biomarkerQC = biomarkerQC.BiomarkerQC(output_folder=test_output_folder, biomarker_outlier_statistics=['iqr'],
                                                   min_number_biomarkers_for_outlier_detection=5,
                                                   min_biomarker_threshold=5,
                                                   max_number_biomarkers_for_outlier_detection=20,
                                                   report_number_of_biomarker_outlier=None, plotting=False)
        self.test_df = pd.DataFrame()
        start = datetime.datetime.today()
        date_list = [start.date() + datetime.timedelta(days=x+7) for x in range(len(np.arange(0, 20, 0.5)) + 1)]
        self.test_df[constant.Columns.DATE.value] = date_list
        values = [np.NAN]
        values.extend(np.arange(0, 20, 0.5))
        # add biomarker flag columns
        for biomarker in constant.Columns.get_biomarker_columns():
            self.test_df[biomarker] = values
            self.test_df[constant.Columns.get_biomarker_flag(biomarker)] = 0
        # add biomarker ratio flag columns
        for biomarker1, biomarker2 in itertools.combinations(constant.Columns.get_biomarker_columns(), 2):
            self.test_df[constant.Columns.get_biomaker_ratio_flag(biomarker1, biomarker2)] = 0

    def tearDown(self) -> None:
        if os.path.exists(test_output_folder):
            shutil.rmtree(test_output_folder)

    def test_biomarker_below_threshold(self):
        flag_biomarker_below_threshold = constant.SewageFlag.BIOMARKER_BELOW_THRESHOLD_OR_EMPTY.value
        correct_flags = []
        correct_flags.extend([flag_biomarker_below_threshold] * 11)  # first three should pass
        correct_flags.extend([0] * 30)  # next five should fail
        self.biomarkerQC.biomarker_below_threshold_or_empty('', self.test_df)
        for biomarker in constant.Columns.get_biomarker_columns():
            valid_result = pd.Series(correct_flags, name=constant.Columns.get_biomarker_flag(biomarker))
            assert_series_equal(self.test_df[constant.Columns.get_biomarker_flag(biomarker)], valid_result)

    def test_calculate_biomarker_ratios(self):
        self.biomarkerQC.biomarker_below_threshold_or_empty('', self.test_df)
        self.biomarkerQC.calculate_biomarker_ratios('', self.test_df)
        valid_ratios = []
        valid_ratios.extend([None for i in range(11)])
        valid_ratios.extend([0.0] * 30)
        for biomarker1, biomarker2 in itertools.combinations(constant.Columns.get_biomarker_columns(), 2):
            valid_result = pd.Series(np.array(valid_ratios, dtype=float), name=biomarker1 + "/" + biomarker2)
            assert_series_equal(self.test_df[biomarker1 + "/" + biomarker2], valid_result)

    def __add_outliers(self):
        idx = 18
        for biomarker in constant.Columns.get_biomarker_columns():
            self.test_df.at[idx, biomarker] = 200  # add outlier
            idx += 2


    def test_detect_outliers(self):
        self.__add_outliers()
        self.biomarkerQC.biomarker_below_threshold_or_empty('', self.test_df)
        self.biomarkerQC.calculate_biomarker_ratios('', self.test_df)
        self.biomarkerQC.detect_outliers('', self.test_df)

        biomarker_ratio_flags = []
        for biomarker1, biomarker2 in itertools.combinations(constant.Columns.get_biomarker_columns(), 2):
            biomarker_ratio_flags.append(constant.Columns.get_biomaker_ratio_flag(biomarker1, biomarker2))
        ## testing ##
        idx = 18
        for biomarker in constant.Columns.get_biomarker_columns():
            biomarker_ratio_flag_columns = [b for b in biomarker_ratio_flags if biomarker in b]
            biomarker_ratio_frame = self.test_df[[b for b in self.test_df.columns if b in biomarker_ratio_flag_columns]]
            biomarker_ratio_flag_row = biomarker_ratio_frame.iloc[idx]
            result = constant.SewageFlag.is_flag_set_for_series(biomarker_ratio_flag_row, constant.SewageFlag.BIOMARKER_RATIO_OUTLIER)
            if not result.all():
                self.fail("Outlier not detected for biomarker: {}".format(biomarker))
            idx += 2

    def test_filter_required_biomarkers(self):
        self.__add_outliers()
        self.biomarkerQC.biomarker_below_threshold_or_empty('', self.test_df)
        self.biomarkerQC.calculate_biomarker_ratios('', self.test_df)
        self.biomarkerQC.detect_outliers('', self.test_df)
       # self.biomarkerQC.analyze_usable_biomarkers('', self.test_df)
        print("here")

# def test_report_last_biomarkers_invalid(self):
#     self.fail()
