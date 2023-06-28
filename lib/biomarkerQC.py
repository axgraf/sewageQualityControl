# Created by alex at 22.06.23
import itertools
from .flags import Flag
from .utils import *
from .plotting import *


class BiomarkerQC:

    def __init__(self, output_folder, outlier_statistics, min_biomarker_threshold,
                 minimal_number_measurements_for_outlier_detection,
                 max_number_measurements_for_outlier_detection, min_number_of_valid_biomarkers_ratios, report_outliers):
        self.output_folder = output_folder
        self.outlier_statistics = outlier_statistics
        self.min_biomarker_threshold = min_biomarker_threshold
        self.minimal_number_measurements_for_outlier_detection = minimal_number_measurements_for_outlier_detection
        self.max_number_measurements_for_outlier_detection = max_number_measurements_for_outlier_detection
        self.min_number_of_valid_biomarkers_ratios = min_number_of_valid_biomarkers_ratios
        self.report_outliers = report_outliers
        self.biomarker_columns = ['biomarker_N1', 'biomarker_N2', 'biomarker_N3', 'biomarker_E', 'biomarker_ORF',
                                  'biomarker_RDRP']
        self.logger = SewageLogger()

    def biomarker_below_threshold(self, sample_location, measurements_df: pd.DataFrame):
        """
        Mark biomarker values that are below given threshold or empty.
        """
        biomarkers_below_dict = dict()
        for biomarker in self.biomarker_columns:
            measurements_df[biomarker + "_below_threshold"] = np.where(
                ((measurements_df[biomarker].isnull()) | (measurements_df[biomarker] < self.min_biomarker_threshold)),
                True, False)
            biomarkers_below_dict[biomarker] = len(
                measurements_df[measurements_df[biomarker + "_below_threshold"] == True])
        self.logger.log.info("[Biomarker below threshold] - [Sample location: '{}'] - \n"
                             "\tNumber of excluded biomarkers which are below threshold ( < {} ) or empty:\n\t\t{}\n"
                             "\tout of {} samples.".
                             format(sample_location, self.min_biomarker_threshold, biomarkers_below_dict,
                                    measurements_df.shape[0]))

    def calculate_biomarker_ratios(self, sample_location, measurements_df: pd.DataFrame):
        """
        Calculate pairwise biomarker ratios only if both are above the minimal threshold level.
        In case one of the biomarker values is below the threshold, the ratio is not calculated.
        """

        biomarkers_ratio_dict = dict()
        for biomarker1, biomarker2 in itertools.combinations(self.biomarker_columns, 2):
            measurements_df[biomarker1 + "/" + biomarker2] = \
                np.where(measurements_df[biomarker1 + "_below_threshold"] | measurements_df[biomarker2 + "_below_threshold"],
                    None,  # if one biomarker is below threshold no threshold will be calculated
                    measurements_df[biomarker1] / measurements_df[biomarker2])
            num_valid_ratios = len(measurements_df[measurements_df[biomarker1 + "/" + biomarker2].notnull()])
            if num_valid_ratios > 0:
                biomarkers_ratio_dict[biomarker1 + "/" + biomarker2] = num_valid_ratios
        self.logger.log.info("[Calculate biomarker ratios] - [Sample location: '{}'] -\n"
                             "\tNumber of calculated biomarker ratios:\n"
                             "\t{}".format(sample_location, biomarkers_ratio_dict))


    def __get_previous_biomarkers_ratios(self, measurements_df: pd.DataFrame, index, biomarker1, biomarker2):
        """
          Obtain last biomarker ratios for the previous measurements starting from current measurement.
          Remove previously detected outliers and empty ratios.
        """
        max_idx = (index - 1 - self.max_number_measurements_for_outlier_detection) \
            if (index - 1 - self.max_number_measurements_for_outlier_detection) >= 0 else 0
        min_idx = index if index >= 0 else 0
        last_N_measurements = measurements_df.iloc[max_idx: min_idx]
        # remove previously detected outliers
        if "outlier_" + biomarker1 + "/" + biomarker2 in last_N_measurements:
            last_N_measurements = last_N_measurements[
                last_N_measurements["outlier_" + biomarker1 + "/" + biomarker2] != True]
        biomarker_ratios = last_N_measurements[
            ['collectionDate', biomarker1 + "/" + biomarker2, biomarker1, biomarker2]]
        # remove empty ratios
        biomarker_ratios = biomarker_ratios.dropna()
        return biomarker_ratios

    def detect_outliers(self, sample_location, measurements_df: pd.DataFrame):
        """
           For each measurement and biomarker pair outliers will be marked.
        """
        self.logger.log.info("[Biomarker outlier detection] - [Sample location: '{}'] - Using following outlier detection methods: {}".format(
                sample_location, self.outlier_statistics))

        stat_dict = dict()
        progress_bar = self.logger.get_progress_bar(measurements_df.shape[0], "Analyzing biomarker outliers:")
        for index, current_measurement in measurements_df.iterrows():
            progress_bar.update(1)
            for biomarker1, biomarker2 in itertools.combinations(self.biomarker_columns, 2):
                add_default_biomarker_statistic(stat_dict, biomarker1, biomarker2)
                biomarker_ratio = biomarker1 + "/" + biomarker2
                if current_measurement[biomarker_ratio]:  # only if biomarker ratio is not None for current measurement
                    stat_dict[biomarker_ratio]["total"] += 1
                    last_biomarker_ratios = self.__get_previous_biomarkers_ratios(measurements_df, index, biomarker1,
                                                                                  biomarker2)
                    measurements_df.at[
                        index, Flag.NOT_ENOUGH_BIOMARKERS_FOR_OUTLIER_DETECTION.name + "_" + biomarker_ratio] = False
                    if len(last_biomarker_ratios) < self.minimal_number_measurements_for_outlier_detection:
                        measurements_df.at[
                            index, Flag.NOT_ENOUGH_BIOMARKERS_FOR_OUTLIER_DETECTION.name + "_" + biomarker_ratio] = True
                        self.logger.log.debug(
                            "[Biomarker outlier detection] - [Sample location: '{}'] - [Collection date: '{}'] - "
                            "Skipping biomarker outlier detection. Found '{}' previous measurements. "
                            "Minimal '{}' previous measurements are required.".format(sample_location,
                                                                                      current_measurement[
                                                                                          'collectionDate'],
                                                                                      len(last_biomarker_ratios),
                                                                                      self.minimal_number_measurements_for_outlier_detection))
                        stat_dict[biomarker1 + "/" + biomarker2]["skipped"] += 1
                        continue
                    is_outlier = self.__is_biomarker_ratio_outlier(current_measurement, last_biomarker_ratios,
                                                                   biomarker1, biomarker2)
                    if is_outlier:
                        stat_dict[biomarker_ratio]["failed"] += 1
                        measurements_df.at[index, "outlier_" + biomarker_ratio] = True
                        measurements_df.at[index, Flag.BIOMARKER_OUTLIER_FLAG.name] = True
                    else:
                        stat_dict[biomarker1 + "/" + biomarker2]["passed"] += 1
                        measurements_df.at[index, "outlier_" + biomarker_ratio] = False
                        measurements_df.at[index, Flag.BIOMARKER_OUTLIER_FLAG.name] = False
                else:
                    self.logger.log.debug(
                        "[Biomarker outlier detection] - [Sample location: '{}'] - [Collection date: '{}'] - "
                        "Empty biomarker ratio '{}' with value: '{}'. Biomarker ratio will be skipped."
                        .format(sample_location, current_measurement['collectionDate'],
                                biomarker_ratio,
                                current_measurement[biomarker_ratio]))

        plot_biomarker_outlier_summary(measurements_df, sample_location,
                                       os.path.join(self.output_folder, "biomarker_qc"), self.outlier_statistics)
        self.logger.log.info("[Biomarker outlier detection] - [Sample location: '{}'] - \n"
                             "Biomarkers ratios:\n{}".format(sample_location,
                                                             pretty_print_biomarker_statistic(stat_dict)))

    def __is_biomarker_ratio_outlier(self, current_measurement, last_biomarker_ratios, biomarker1, biomarker2):
        """
            Detect outliers for a given biomarker ratio based on previous measurements
            using random forest and/or interquartile range and/or z-score outlier detection methods
        """
        biomarker_ratio = biomarker1 + "/" + biomarker2
        current_biomarker_ratio = current_measurement[biomarker_ratio]
        outlier_detected = []
        if 'lof' in self.outlier_statistics or 'all' in self.outlier_statistics:
            is_lof_outlier = is_outlier_local_outlier_factor(current_biomarker_ratio, last_biomarker_ratios[biomarker_ratio])
            outlier_detected.append(is_lof_outlier)
        if 'rf' in self.outlier_statistics or 'all' in self.outlier_statistics:
            is_isolation_forest_outlier, isolation_forest_score_ratio = is_outlier_isolation_forest(
                current_biomarker_ratio, last_biomarker_ratios[biomarker_ratio])
            outlier_detected.append(is_isolation_forest_outlier)
        #        is_confidence_interval_99_outlier, confidence_interval = is_confidence_interval_outlier(current_biomarker_ratio,
        #                last_biomarker_ratios[biomarker_ratio],  0.99)
        if 'iqr' in self.outlier_statistics or 'all' in self.outlier_statistics:
            is_iqr_outlier, iqr_range = inter_quantil_range(current_biomarker_ratio,
                                                            last_biomarker_ratios[biomarker_ratio])
            outlier_detected.append(is_iqr_outlier)
        if 'zscore' in self.outlier_statistics or 'all' in self.outlier_statistics:
            is_zscore_outlier, z_score_threshold = is_outlier_modified_z_score(current_biomarker_ratio,
                                                                               last_biomarker_ratios[biomarker_ratio])
            outlier_detected.append(is_zscore_outlier)
        if all(outlier_detected):
            return True
        return False

    def filter_required_biomarkers(self, sample_location, measurements_df: pd.DataFrame):
        """
        Sums the number of biomarkers which are not marked as outliers using the outlier columns.
        Flags measurments which have less than specified number of valid biomarkers.
        """
        biomarker_columns = []
        for biomarker1, biomarker2 in itertools.combinations(self.biomarker_columns, 2):
            biomarker_columns.append(biomarker1 + "/" + biomarker2)
        biomarker_frame = measurements_df.loc[:, [i in biomarker_columns for i in measurements_df.columns]]
        measurements_df['biomarkerRatioSums'] = biomarker_frame.count(axis=1)
        measurements_df[Flag.MIN_BIOMARKER_NUMBER_NOT_REACHED.name] = \
            np.where(measurements_df['biomarkerRatioSums'] >= self.min_number_of_valid_biomarkers_ratios, False, True)

        count_biomarker_ready_for_analysis = measurements_df[Flag.MIN_BIOMARKER_NUMBER_NOT_REACHED.name].value_counts()
        self.logger.log.info("[Select Biomarkers] - [Sample location: '{}'] - \n"
                             "Samples with >= '{}' biomarker ratios:\n{}".format(sample_location,
                                                                                 self.min_number_of_valid_biomarkers_ratios,
                                                                                 count_biomarker_ready_for_analysis.to_string()))

    def report_last_biomarkers_invalid(self, sample_location, measurements_df: pd.DataFrame):
        last_two_measurements = measurements_df.tail(2)
        for biomarker1, biomarker2 in itertools.combinations(self.biomarker_columns, 2):
            outlier_column_name = "outlier_" + biomarker1 + "/" + biomarker2
            if outlier_column_name in last_two_measurements.columns:
                false_biomarkers = last_two_measurements[
                    last_two_measurements["outlier_" + biomarker1 + "/" + biomarker2] == True]
                number_of_outliers = false_biomarkers.shape[0]
                if number_of_outliers >= self.report_outliers:
                    biomarker_last_outliers = last_two_measurements[
                        ['collectionDate', biomarker1, biomarker2, biomarker1 + "/" + biomarker2]]
                    self.logger.log.warn("[Report last two biomarker outliers] - [Sample location: '{}'] - "
                                         "The last two measurements for biomarker ratio '{}/{}' are marked as outliers.\n"
                                         "{}".format(sample_location, biomarker1, biomarker2, biomarker_last_outliers))
                    # Todo: handle reporting
