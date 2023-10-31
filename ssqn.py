#!/usr/bin/env python3

import os
import itertools
import argparse
import pickle

import numpy as np
import pandas as pd

from lib.arcgis import *
from matplotlib.backends.backend_pdf import PdfPages
from lib.config import Config
from lib.constant import *
from lib.biomarkerQC import BiomarkerQC
from lib.surrogatevirusQC import SurrogateVirusQC
from lib.water_quality import WaterQuality
from lib.sewage_flow import SewageFlow
from lib.normalization import SewageNormalization
import lib.arcgis as arcgis
import lib.bayvoc as bayVocConn
import lib.utils as utils
import lib.statistics as sewageStat
import lib.plotting as plotting
import lib.database as db
import re


class SewageQuality:

    def __init__(self, input_file, config_file, output_folder, verbosity, quiet, rerun_all, biomarker_outlier_statistics,
                 min_biomarker_threshold, min_number_biomarkers_for_outlier_detection,
                 max_number_biomarkers_for_outlier_detection, report_number_of_biomarker_outlier, periode_month_surrogatevirus,
                 surrogatevirus_outlier_statistics, min_number_surrogatevirus_for_outlier_detection,
                 water_quality_number_of_last_month, min_number_of_last_measurements_for_water_qc, water_qc_outlier_statistics,
                 fraction_last_samples_for_dry_flow, min_num_samples_for_mean_dry_flow, heavy_precipitation_factor,
                 mean_sewage_flow_below_typo_factor, mean_sewage_flow_above_typo_factor, min_number_of_biomarkers_for_normalization,
                 base_reproduction_value_factor, num_previous_days_reproduction_factor, max_number_of_flags_for_outlier
                  ):

        self.input_file = input_file
        self.config_file = config_file
        self.sewage_samples = None
        self.output_folder = output_folder
        self.verbosity = verbosity
        self.quiet = quiet
        self.rerun_all = rerun_all
        # biomarker qc
        self.biomarker_outlier_statistics = biomarker_outlier_statistics
        self.min_biomarker_threshold = min_biomarker_threshold
        self.min_number_biomarkers_for_outlier_detection = min_number_biomarkers_for_outlier_detection
        self.max_number_biomarkers_for_outlier_detection = max_number_biomarkers_for_outlier_detection
        self.report_number_of_biomarker_outlier = report_number_of_biomarker_outlier
        # Surrogate virus
        self.periode_month_surrogatevirus = periode_month_surrogatevirus
        # Sewage flow
        self.fraction_last_samples_for_dry_flow = fraction_last_samples_for_dry_flow
        self.min_num_samples_for_mean_dry_flow = min_num_samples_for_mean_dry_flow
        self.heavy_precipitation_factor = heavy_precipitation_factor
        self.mean_sewage_flow_below_typo_factor = mean_sewage_flow_below_typo_factor
        self.mean_sewage_flow_above_typo_factor = mean_sewage_flow_above_typo_factor
        self.surrogatevirus_outlier_statistics = surrogatevirus_outlier_statistics
        self.min_number_surrogatevirus_for_outlier_detection = min_number_surrogatevirus_for_outlier_detection
        # Water quality
        self.water_quality_number_of_last_month = water_quality_number_of_last_month
        self.min_number_of_last_measurements_for_water_qc = min_number_of_last_measurements_for_water_qc
        self.water_qc_outlier_statistics = water_qc_outlier_statistics
        # biomarker normalization
        self.min_number_of_biomarkers_for_normalization = min_number_of_biomarkers_for_normalization
        self.base_reproduction_value_factor = base_reproduction_value_factor
        self.num_previous_days_reproduction_factor = num_previous_days_reproduction_factor
        self.max_number_of_flags_for_outlier = max_number_of_flags_for_outlier
        self.sewageStat = sewageStat.SewageStat()
        self.logger = utils.SewageLogger(self.output_folder, verbosity=verbosity, quiet=quiet)
        self.__load_data()
        self.__initialize()

    def __load_data(self):
        if self.input_file:
            self.sewage_samples_dict = utils.read_excel_input_files(self.input_file)
        elif self.config_file:
            config = Config(self.config_file)
#            arcgis = Arcgis(config)
#            sewage_samples_dict = arcgis.obtain_sewage_samples()
#            with open("data/sewageData_arcgis.dat", "wb") as writer:
#                pickle.dump(sewage_samples_dict, writer)
            sewage_samples_dict = None
            with open('data/sewageData_arcgis.dat', 'rb') as f:
                sewage_samples_dict = pickle.load(f)
            self.sewage_samples_dict = dict()
            for loc, samples_list in sewage_samples_dict.items():
                if 'A-STADT' in loc:
                    self.sewage_samples_dict[loc] = pd.DataFrame(s.__dict__ for s in samples_list)
           # bayvoc = bayVocConn.BayVOC(config)
           # self.sewage_samples_dict = bayvoc.read_all_sewagesamples_from_db()
        # arcgis = Arcgis(Config(self.config))
        # self.sewage_samples = arcgis.obtain_sewage_samples()
        # Todo: switch to real data import
#        import pickle
#        with open('data/sewageData.dat', 'rb') as f:
#            self.sewage_samples = pickle.load(f)
#        with open('data/sewagePlantData.dat', 'rb') as f:
#            self.sewage_plants2trockenwetterabfluss = pickle.load(f)
        self.sewage_plants2trockenwetterabfluss = dict()

    def __initialize(self):
        if not os.path.exists(self.output_folder):
            os.makedirs(self.output_folder)
        self.database = db.SewageDatabase()
        self.biomarkerQC = BiomarkerQC(self.output_folder, self.sewageStat, self.biomarker_outlier_statistics, self.min_biomarker_threshold,
                                  self.min_number_biomarkers_for_outlier_detection,
                                  self.max_number_biomarkers_for_outlier_detection,
                                  self.report_number_of_biomarker_outlier)
        self.water_quality = WaterQuality(self.output_folder, self.sewageStat, self.water_quality_number_of_last_month,
                                          self.min_number_of_last_measurements_for_water_qc, self.water_qc_outlier_statistics)
        self.sewage_flow = SewageFlow(self.output_folder, self.sewageStat, self.sewage_plants2trockenwetterabfluss,
                                      self.fraction_last_samples_for_dry_flow, self.min_num_samples_for_mean_dry_flow,
                                      self.heavy_precipitation_factor, self.mean_sewage_flow_below_typo_factor, self.mean_sewage_flow_above_typo_factor)
        self.surrogateQC = SurrogateVirusQC(self.sewageStat, self.periode_month_surrogatevirus,
                                     self.min_number_surrogatevirus_for_outlier_detection,
                                     self.biomarker_outlier_statistics, self.output_folder)
        self.sewageNormalization = SewageNormalization(self.sewageStat, self.max_number_of_flags_for_outlier, self.min_number_of_biomarkers_for_normalization,
                                                       self.base_reproduction_value_factor, self.num_previous_days_reproduction_factor, self.output_folder)



    def __initalize_columns(self, measurements: pd.DataFrame):
        for biomarker in Columns.get_biomarker_columns():
            measurements[CalculatedColumns.get_biomarker_flag(biomarker)] = 0
        for biomarker1, biomarker2 in itertools.combinations(Columns.get_biomarker_columns(), 2):
            measurements[biomarker1 + "/" + biomarker2] = np.NAN
            measurements[CalculatedColumns.get_biomaker_ratio_flag(biomarker1, biomarker2)] = 0
        for c in CalculatedColumns:
            if not c.value in measurements:
                if c.type == bool:
                    measurements[c.value] = False
                elif c.type == str:
                    measurements[c.value] = ""
                else:
                    measurements[c.value] = 0
                measurements[c.value] = measurements[c.value].astype(c.type)

    def save_dataframe(self, sample_location, measurements: pd.DataFrame):
        result_folder = os.path.join(self.output_folder, "results")
        if not os.path.exists(result_folder):
            os.makedirs(result_folder)
        output_file = os.path.join(result_folder, "normalized_sewage_{}.xlsx".format(sample_location))
        measurements.to_excel(output_file, index=False)

    def __setup(self, sample_location, measurements: pd.DataFrame):
        plausibility_dict = {}
        measurements = measurements.fillna(value=np.nan)
        measurements[Columns.DATE.value] = measurements[Columns.DATE.value].replace({r'(\d+-\d+-\d+).*': r'\1'}, regex=True)

        # check if the date format is YYYY-mm-dd
        format_not_plausible = measurements[measurements[Columns.DATE.value].str.contains("[\d{4}\-\d{2}\-\d{2}]").eq(False)].index.tolist()
        if len(format_not_plausible) > 0:
            plausibility_dict["format_not_plausible_index"] = format_not_plausible
        measurements[Columns.DATE.value] = pd.to_datetime(measurements[Columns.DATE.value], format="%Y-%m-%d").dt.normalize()

        # check if the year of the measurement is plausible: Corona pandemic started in 2020
        year_not_plausible = measurements[(measurements[Columns.DATE.value].dt.year > pd.to_datetime("2019", format="%Y").year).eq(False)].index.tolist()
        if len(year_not_plausible) > 0:
            plausibility_dict["year_not_plausible_index"] = year_not_plausible

        # check if the month value is <=12
        month_not_plausible = measurements[(measurements[Columns.DATE.value].dt.month <= pd.to_datetime("12", format="%m").month).eq(False)].index.tolist()
        if len(month_not_plausible) > 0:
            plausibility_dict["month_not_plausible_index"] = month_not_plausible

        # check if the month value is <=31
        day_not_plausible = measurements[(measurements[Columns.DATE.value].dt.day <= pd.to_datetime("31", format="%d").day).eq(False)].index.tolist()
        if len(day_not_plausible) > 0:
            plausibility_dict["day_not_plausible_index"] = day_not_plausible

        # setup comment fields
        measurements[[Columns.COMMENT_ANALYSIS.value, Columns.COMMENT_OPERATION.value]] = \
            measurements[[Columns.COMMENT_ANALYSIS.value, Columns.COMMENT_OPERATION.value]].apply(
                lambda x: x.astype(str).str.strip()
            )
        measurements[[Columns.COMMENT_ANALYSIS.value, Columns.COMMENT_OPERATION.value]] = \
            measurements[[Columns.COMMENT_ANALYSIS.value, Columns.COMMENT_OPERATION.value]].replace('nan', '')
        measurements[[Columns.TROCKENTAG.value]] = measurements[[Columns.TROCKENTAG.value]].astype(str)
        # measurements = measurements.drop(columns=["flags"])   # artefact from stored data --> will be removed later
        # measurements = utils.convert_sample_list2pandas(measurements)
        # Sort by collection date. Newest last.
        measurements.sort_values(by=Columns.DATE.value, ascending=True, inplace=True, ignore_index=True)
        self.__initalize_columns(measurements)
        self.database.needs_recalcuation(sample_location, measurements, self.rerun_all)
        return plausibility_dict, measurements

    def __is_plot_not_generated(self, sample_location):
        return not os.path.exists(os.path.join(self.output_folder, "plots", "{}.plots.pdf".format(sample_location)))

    def __plot_results(self, measurements: pd.DataFrame, sample_location):
        if not os.path.exists(os.path.join(self.output_folder, "plots")):
            os.makedirs(os.path.join(self.output_folder, "plots"))
        pdf_pages = PdfPages(os.path.join(self.output_folder, "plots", "{}.plots.pdf".format(sample_location)))
        plotting.plot_biomarker_outlier_summary(pdf_pages, measurements, sample_location, self.biomarker_outlier_statistics)
        plotting.plot_surrogatvirus(pdf_pages, measurements, sample_location, self.surrogatevirus_outlier_statistics)
        plotting.plot_sewage_flow(pdf_pages, measurements, sample_location)
        plotting.plot_water_quality(pdf_pages, measurements, sample_location, self.water_qc_outlier_statistics)
        plotting.plot_biomarker_normalization(pdf_pages, measurements, sample_location)
        plotting.plot_general_outliers(pdf_pages, measurements, sample_location)
        pdf_pages.close()

    def run_quality_control(self):
        """
        Main method to run the quality checks and normalization
        """
        for idx, (sample_location, measurements) in enumerate(self.sewage_samples_dict.items()):
            self.logger.log.info("\n####################################################\n"
                                 "\tSewage location: {} "
                                 "\n####################################################".format(sample_location))
            plausibility_dict, measurements = self.__setup(sample_location, measurements)
            ### Plausibilitätscheck: dict with the index of the odd values
            if len(plausibility_dict) > 0:
                self.logger.log.info("Check date filed:{}".format(plausibility_dict))
            self.logger.log.info("{}/{} new measurements to analyze".format(CalculatedColumns.get_num_of_unprocessed(measurements),
                                                                            measurements.shape[0]))
            progress_bar = self.logger.get_progress_bar(CalculatedColumns.get_num_of_unprocessed(measurements), "Analyzing samples")
            self.sewageStat.set_sample_location_and_total_number(sample_location, CalculatedColumns.get_num_of_unprocessed(measurements))
            changes_detected = False
            for index, current_measurement in measurements.iterrows():
                if CalculatedColumns.needs_processing(current_measurement):
                    changes_detected = True
                    progress_bar.update(1)
                    # -----------------  BIOMARKER QC -----------------------
                    # 1. check for comments. Flag samples that contain any commentary.
                    self.biomarkerQC.check_comments(sample_location, measurements, index)
                    self.biomarkerQC.check_mean_sewage_flow_present(sample_location, measurements, index)
                    # 2. Mark biomarker values below threshold which are excluded from the analysis.
                    self.biomarkerQC.biomarker_below_threshold_or_empty(sample_location, measurements, index)
                    # 3. Calculate pairwise biomarker values if biomarkers were not marked to be below threshold.
                    self.biomarkerQC.calculate_biomarker_ratios(sample_location, measurements, index)
                    # 4. Detect outliers
                    self.biomarkerQC.detect_outliers(sample_location, measurements, index)
                    # 5. Assign biomarker outliers based on ratio outliers
                    self.biomarkerQC.assign_biomarker_outliers_based_on_ratio_flags(sample_location, measurements, index)
                    self.biomarkerQC.analyze_usable_biomarkers(sample_location, measurements, index)
                    # 6. Create report in case the last two biomarkers were identified as outliers
                    # self.biomarkerQC.report_last_biomarkers_invalid(sample_location, measurements)

                    # --------------------  SUROGATVIRUS QC -------------------
                    self.surrogateQC.filter_dry_days_time_frame(sample_location, measurements, index)
                    self.surrogateQC.is_surrogatevirus_outlier(sample_location, measurements, index)

                    # --------------------  SEWAGE FLOW -------------------
                    self.sewage_flow.sewage_flow_quality_control(sample_location, measurements, index)

                    # --------------------  WATER QUALITY -------------------
                    self.water_quality.check_water_quality(sample_location, measurements, index)

                    # --------------------  NORMALIZATION -------------------
                    self.sewageNormalization.normalize_biomarker_values(sample_location, measurements, index)

                    # --------------------  MARK OUTLIERS FROM ALL STEPS -------------------
                    self.sewageNormalization.decide_biomarker_usable_based_on_flags(sample_location, measurements, index)
            progress_bar.close()
            print("    ")
            if changes_detected or self.__is_plot_not_generated(sample_location):
                self.logger.log.info(self.sewageStat.print_statistics())
                self.logger.log.info("Generating plots...")
                self.__plot_results(measurements, sample_location)
                # Experimental: Final step explain flags
                measurements['flags_explained'] = SewageFlag.explain_flag_series(measurements[CalculatedColumns.FLAG.value])
                self.logger.log.info("Add '{}' to database...".format(sample_location))
                self.database.add_sewage_location2db(sample_location, measurements)
                self.logger.log.info("Export '{}' to excel file...".format(sample_location))
                self.save_dataframe(sample_location, measurements)



if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Sewage qPCR quality control",
        usage='use "python3 ssqn.py --help" for more information',
        epilog="author: Dr. Alexander Graf (graf@genzentrum.lmu.de)", formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-i', '--input', metavar="FILE", type=str,
                        help="Specifiy input excel file with biomarker values",
                        required=False)
    parser.add_argument('-c', '--config', metavar="FILE", type=str,
                        help="Config file for DB connections",
                        required=False)
    parser.add_argument('-o', '--output_folder', metavar="FOLDER", default="sewage_qc", type=str,
                        help="Specifiy output folder. (default folder: 'sewage_qc')",
                        required=False)
    parser.add_argument('-r', '--rerun_all', action="store_true", help="Rerun the analysis on all samples.")
    parser.add_argument('-v', '--verbosity', action="count", help="Increase output verbosity.")
    parser.add_argument('-q', '--quiet', action='store_true', help="Print litte output.")

    biomarker_qc_group = parser.add_argument_group("Biomarker quality control")
    biomarker_qc_group.add_argument('--biomarker_outlier_statistics', metavar="METHOD", default=['iqr', 'lof'], nargs='+',
                        help=("Which outlier detection methods should be used? Multiple selections allowed. (default: 'lof','iqr')\n"
                              "Possible choices are : [lof, rf, iqr, zscore, ci, all]\n"
                              "E.g. to select 'rf' and 'iqr' use: --biomarker_outlier_statistics rf iqr \n"
                              "\tsvm = one class SVM\n"
                              "\tlof = local outlier factor\n"
                              "\trf = random forest\n"
                              "\tiqr = interquartile range\n"
                              "\tzscore = modified z-score\n"
                              "\tci = 99%% confidence interval\n"
                              "\tall = use all methods\n") ,
                        choices=["lof", "rf", "iqr", "zscore", "ci", "svm", "all"],
                        required=False)
    biomarker_qc_group.add_argument('--biomarker_min_threshold', metavar="FLOAT", default=1.5, type=float,
                        help="Minimal biomarker threshold. (default: 1.5)",
                        required=False)
    biomarker_qc_group.add_argument('--min_number_biomarkers_for_outlier_detection', metavar="INT", default=9, type=int,
                        help="Minimal number of previous measurements required for outlier detection, otherwise this step is skipped. (default: 9)",
                        required=False)
    biomarker_qc_group.add_argument('--max_number_biomarkers_for_outlier_detection', metavar="INT", default=50, type=int,
                        help="Maximal number of previous measurements to use for outlier detection. (default: 50)",
                        required=False)
    biomarker_qc_group.add_argument('--report_number_of_biomarker_outliers', metavar="INT", default=2, type=int,
                        help="The number of outliers identified in the last N consecutive biomarker ratios that trigger a report. (default: 2)",
                        required=False)

    surrogatevirus_group = parser.add_argument_group("Surrogate virus quality control")
    surrogatevirus_group.add_argument('--periode_month_surrogatevirus', metavar="INT", default=4, type=int,
                        help="The periode of time (month) taken into account for surrogatevirus outliers. (default: 4)",
                        required=False)
    surrogatevirus_group.add_argument('--min_number_surrogatevirus_for_outlier_detection', metavar="INT", default=9, type=int,
                        help="Minimal number of surrogatevirus measurements. (default: 9)",
                        required=False)
    surrogatevirus_group.add_argument('--surrogatevirus_outlier_statistics', metavar="METHOD", default=['iqr', 'lof'], nargs='+',
                                    help=(
                                        "Which outlier detection methods should be used for surrogatevirus qc? Multiple selections allowed. (default: 'lof','iqr')\n"
                                        "Possible choices are : [lof, rf, iqr, zscore, ci, all]\n"
                                        "E.g. to select 'rf' and 'iqr' use: --outlier_statistics rf iqr \n"
                                        "\tsvm = one class SVM\n"
                                        "\tlof = local outlier factor\n"
                                        "\trf = random forest\n"
                                        "\tiqr = interquartile range\n"
                                        "\tzscore = modified z-score\n"
                                        "\tci = 99%% confidence interval\n"
                                        "\tall = use all methods\n"),
                                    choices=["lof", "rf", "iqr", "zscore", "ci", "svm", "all"],
                                    required=False)

    sewage_flow_group = parser.add_argument_group("Sewage flow quality control")
    sewage_flow_group.add_argument('--fraction_last_samples_for_dry_flow', metavar="FLOAT", default=0.1, type=float,
                                     help="If the dry flow of the sewage treatment plant is not known, "
                                          "the average dry flow rate is estimated from the previous flows rate of "
                                          "the lowest N percent of the samples. (default: 0.1)", choices=[round(x * 0.1, 1) for x in range(0, 10)],
                                     required=False)
    sewage_flow_group.add_argument('--min_num_samples_for_mean_dry_flow', metavar="INT", default=5, type=int,
                                   help="If the dry flow of the sewage treatment plant is not known, "
                                        "minimal N previous samples are required for the estimation of the dry flow rate. (default: 5)",
                                   required=False)
    sewage_flow_group.add_argument('--heavy_precipitation_factor', metavar="FLOAT", default=2.0, type=float,
                                   help="Factor above which the mean flow must be in comparison to the dry weather "
                                        "flow in order for the sample to be sorted out as a heavy rain event. (default: 2.0)",
                                   required=False)
    sewage_flow_group.add_argument('--mean_sewage_flow_below_typo_factor', metavar="FLOAT", default=1.5, type=float,
                                   help="Factor below which the mean flow must be in comparison to the dry weather "
                                        "flow in order to mark the value as a probable typo. (default: 1.5)",
                                   required=False)
    sewage_flow_group.add_argument('--mean_sewage_flow_above_typo_factor', metavar="FLOAT", default=9.0, type=float,
                                   help="Factor above which the mean flow must be in comparison to the dry weather "
                                        "flow in order to mark the value as a probable typo. (default: 9.0)",
                                   required=False)

    water_quality_group = parser.add_argument_group("Water quality control")
    water_quality_group.add_argument('--water_quality_number_of_last_month', metavar="INT", default=4, type=int,
                                      help="The number of last months to be used for water quality testing. (default: 4)",
                                      required=False)
    water_quality_group.add_argument('--min_number_of_last_measurements_for_water_qc', metavar="INT", default=9, type=int,
                                     help="The minimal number of last measurements required for water quality quality control. (default: 9)",
                                     required=False)
    water_quality_group.add_argument('--water_qc_outlier_statistics', metavar="METHOD", default=['iqr', 'lof'], nargs='+',
                                    help=(
                                        "Which outlier detection methods should be used for water qc? Multiple selections allowed. (default: 'lof','rf','iqr')\n"
                                        "Possible choices are : [lof, rf, iqr, zscore, ci, all]\n"
                                        "E.g. to select 'rf' and 'iqr' use: --outlier_statistics rf iqr \n"
                                        "\tsvm = one class SVM\n"
                                        "\tlof = local outlier factor\n"
                                        "\trf = random forest\n"
                                        "\tiqr = interquartile range\n"
                                        "\tzscore = modified z-score\n"
                                        "\tci = 99%% confidence interval\n"
                                        "\tall = use all methods\n"),
                                    choices=["lof", "rf", "iqr", "zscore", "ci", "svm", "all"],
                                    required=False)
    reproduction_group = parser.add_argument_group("Reproduction factor")
    reproduction_group.add_argument('--num_previous_days_reproduction_factor', metavar="INT", default=7, type=int,
                                    help="Number of days to include samples for outlier analysis of the reproduction factor. (default: 7)",
                                    required=False)
    reproduction_group.add_argument('--base_reproduction_value_factor', metavar="FLOAT", default=3.8, type=float,
                                     help="Factor below which the mean normalized biomarker value must be in comparison "
                                          "to the mean normalized biomarker value of the last 7 days. (default: 3.8)",
                                     required=False)
    normalization_group = parser.add_argument_group("Biomarker normalization")
    normalization_group.add_argument('--max_number_of_flags_for_outlier', metavar="INT", default=2, type=int,
                                     help="Maximal number of accumulated flags from all quality controls. "
                                          "If the number is higher the sample will be marked as outlier. (default: 2)",
                                     required=False)
    normalization_group.add_argument('--min_number_of_biomarkers_for_normalization', metavar="INT", default=2, type=int,
                                     help="Minimal number of biomarkers used for normalization. (default: 2)",
                                     required=False)


    args = parser.parse_args()
    sewageQuality = SewageQuality(args.input, args.config, args.output_folder, args.verbosity, args.quiet, args.rerun_all,
                                  args.biomarker_outlier_statistics, args.biomarker_min_threshold,
                                  args.min_number_biomarkers_for_outlier_detection,
                                  args.max_number_biomarkers_for_outlier_detection,
                                  args.report_number_of_biomarker_outliers, args.periode_month_surrogatevirus,
                                  args.surrogatevirus_outlier_statistics, args.min_number_surrogatevirus_for_outlier_detection,
                                  args.water_quality_number_of_last_month,
                                  args.min_number_of_last_measurements_for_water_qc, args.water_qc_outlier_statistics,
                                  args.fraction_last_samples_for_dry_flow, args.min_num_samples_for_mean_dry_flow,
                                  args.heavy_precipitation_factor, args.mean_sewage_flow_below_typo_factor,
                                  args.mean_sewage_flow_above_typo_factor, args.min_number_of_biomarkers_for_normalization,
                                  args.base_reproduction_value_factor, args.num_previous_days_reproduction_factor, args.max_number_of_flags_for_outlier)

    sewageQuality.run_quality_control()

