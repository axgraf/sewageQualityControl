# Created by alex at 22.06.23
import math
from .utils import *
from .plotting import *


class BiomarkerQC:

    def __init__(self, output_folder, biomarker_outlier_statistics, min_biomarker_threshold,
                 min_number_biomarkers_for_outlier_detection,
                 max_number_biomarkers_for_outlier_detection, report_number_of_biomarker_outlier):
        self.output_folder = output_folder
        self.biomarker_outlier_statistics = biomarker_outlier_statistics
        self.min_biomarker_threshold = min_biomarker_threshold
        self.min_number_biomarkers_for_outlier_detection = min_number_biomarkers_for_outlier_detection
        self.max_number_biomarkers_for_outlier_detection = max_number_biomarkers_for_outlier_detection
        self.report_number_of_biomarker_outlier = report_number_of_biomarker_outlier
        self.logger = SewageLogger()

    def biomarker_below_threshold(self, sample_location, measurements_df: pd.DataFrame):
        """
        Mark biomarker values that are below given threshold or empty.
        """
        biomarkers_below_dict = dict()
        flag_value = SewageFlag.BIOMARKER_BELOW_THRESHOLD.value
        for biomarker in Columns.get_biomarker_columns():
            measurements_df[Columns.get_biomarker_flag(biomarker)] += np.where(
                ((measurements_df[biomarker].isnull()) | (measurements_df[biomarker] < self.min_biomarker_threshold)),
                flag_value, 0)
            biomarker_flag_series = measurements_df[Columns.get_biomarker_flag(biomarker)]
            biomarkers_below_dict[biomarker] = \
                len(measurements_df[SewageFlag.is_flag_set_for_series(biomarker_flag_series, SewageFlag.BIOMARKER_BELOW_THRESHOLD)])
        self.logger.log.info("[Biomarker below threshold] - [Sample location: '{}'] - \n"
                             "\tNumber of excluded biomarkers which are below threshold ( < {} ) or empty:\n\t{}\n"
                             "\tout of {} samples.".
                             format(sample_location, self.min_biomarker_threshold, biomarkers_below_dict,
                                    measurements_df.shape[0]))

    def calculate_biomarker_ratios(self, sample_location, measurements_df: pd.DataFrame):
        """
        Calculate pairwise biomarker ratios only if both are above the minimal threshold level.
        In case one of the biomarker values is below the threshold, the ratio is not calculated.
        """
        biomarkers_ratio_dict = dict()
        for biomarker1, biomarker2 in itertools.combinations(Columns.get_biomarker_columns(), 2):
            biomarker1_series = measurements_df[Columns.get_biomarker_flag(biomarker1)]
            biomarker2_series = measurements_df[Columns.get_biomarker_flag(biomarker2)]
            measurements_df[biomarker1 + "/" + biomarker2] = \
                np.where((SewageFlag.is_flag_set_for_series(biomarker1_series, SewageFlag.BIOMARKER_BELOW_THRESHOLD)) |
                         (SewageFlag.is_flag_set_for_series(biomarker2_series, SewageFlag.BIOMARKER_BELOW_THRESHOLD)),
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
        max_idx = (index - 1 - self.max_number_biomarkers_for_outlier_detection) \
            if (index - 1 - self.max_number_biomarkers_for_outlier_detection) >= 0 else 0
        min_idx = index if index >= 0 else 0
        last_N_measurements = measurements_df.iloc[max_idx: min_idx]
        # remove previously detected outliers
        biomarker_ratio_flags = last_N_measurements[Columns.get_biomaker_ratio_flag(biomarker1, biomarker2)]
        last_N_measurements = last_N_measurements[SewageFlag.is_not_flag_set_for_series(biomarker_ratio_flags,
                                                                                        SewageFlag.BIOMARKER_RATIO_OUTLIER)]
        biomarker_ratios = last_N_measurements[[Columns.DATE.value, biomarker1 + "/" + biomarker2, biomarker1, biomarker2]]
        # remove empty ratios
        biomarker_ratios = biomarker_ratios.dropna()
        return biomarker_ratios

    def detect_outliers(self, sample_location, measurements_df: pd.DataFrame):
        """
           For each measurement and biomarker pair outliers will be marked.
        """
        self.logger.log.info("[Biomarker outlier detection] - [Sample location: '{}'] - Using following outlier detection methods: {}".format(
                sample_location, self.biomarker_outlier_statistics))
        stat_dict = dict()
        progress_bar = self.logger.get_progress_bar(measurements_df.shape[0], "Analyzing biomarker outliers")
        for index, current_measurement in measurements_df.iterrows():
            progress_bar.update(1)
            for biomarker1, biomarker2 in itertools.combinations(Columns.get_biomarker_columns(), 2):
                add_default_biomarker_statistic(stat_dict, biomarker1, biomarker2)
                biomarker_ratio = biomarker1 + "/" + biomarker2
                if current_measurement[biomarker_ratio] and not math.isnan(current_measurement[biomarker_ratio]):  # only if biomarker ratio is not None for current measurement
                    stat_dict[biomarker_ratio]["total"] += 1
                    last_biomarker_ratios = self.__get_previous_biomarkers_ratios(measurements_df, index, biomarker1,
                                                                                  biomarker2)
                    if len(last_biomarker_ratios) < self.min_number_biomarkers_for_outlier_detection:
                        measurements_df.at[
                            index, Columns.get_biomaker_ratio_flag(biomarker1, biomarker2)] += SewageFlag.NOT_ENOUGH_PREVIOUS_BIOMARKER_VALUES.value
                        self.logger.log.debug(
                            "[Biomarker outlier detection] - [Sample location: '{}'] - [Collection date: '{}'] - "
                            "Skipping biomarker outlier detection. Found '{}' previous measurements. "
                            "Minimal '{}' previous measurements are required.".format(sample_location, current_measurement[Columns.DATE.value],
                                                                                      len(last_biomarker_ratios),
                                                                                      self.min_number_biomarkers_for_outlier_detection))
                        stat_dict[biomarker1 + "/" + biomarker2]["skipped"] += 1
                        continue
                    is_outlier = detect_outliers(self.biomarker_outlier_statistics, last_biomarker_ratios[biomarker_ratio], current_measurement[biomarker_ratio])
                    if is_outlier:
                        stat_dict[biomarker_ratio]["failed"] += 1
                        measurements_df.at[index, Columns.get_biomaker_ratio_flag(biomarker1, biomarker2)] += SewageFlag.BIOMARKER_RATIO_OUTLIER.value
                        measurements_df.at[index, Columns.get_biomarker_flag(biomarker1)] += SewageFlag.BIOMARKER_OUTLIER.value
                        measurements_df.at[index, Columns.get_biomarker_flag(biomarker2)] += SewageFlag.BIOMARKER_OUTLIER.value
                    else:
                        stat_dict[biomarker1 + "/" + biomarker2]["passed"] += 1
                else:
                    stat_dict[biomarker1 + "/" + biomarker2]["skipped"] += 1
                    self.logger.log.debug(
                        "[Biomarker outlier detection] - [Sample location: '{}'] - [Collection date: '{}'] - "
                        "Empty biomarker ratio '{}' with value: '{}'. Biomarker ratio will be skipped."
                        .format(sample_location, current_measurement[Columns.DATE.value],
                                biomarker_ratio,
                                current_measurement[biomarker_ratio]))

        plot_biomarker_outlier_summary(measurements_df, sample_location,
                                       os.path.join(self.output_folder, "biomarker"), self.biomarker_outlier_statistics)
        self.logger.log.info("[Biomarker outlier detection] - [Sample location: '{}'] - \n"
                             "Biomarkers ratios:\n{}".format(sample_location,
                                                             pretty_print_biomarker_statistic(stat_dict)))


    def filter_required_biomarkers(self, sample_location, measurements_df: pd.DataFrame):
        """
        Sums the number of biomarkers which are not marked as outliers using the outlier columns.
        Flags measurments which have less than specified number of valid biomarkers.
        """
        biomarker_flags_columns = []
        for biomarker in Columns.get_biomarker_columns():
            biomarker_flags_columns.append(Columns.get_biomarker_flag(biomarker))
        biomarker_flags_frame = measurements_df[[c for c in measurements_df.columns if c in biomarker_flags_columns]]
        valid_biomarkers = SewageFlag.is_not_flag_set_for_series(biomarker_flags_frame, SewageFlag.BIOMARKER_OUTLIER | SewageFlag.BIOMARKER_BELOW_THRESHOLD)
        measurements_df[Columns.NUMBER_OF_USABLE_BIOMARKERS.value] = valid_biomarkers.sum(axis=1)
        flag_biomarker_outlier = SewageFlag.BIOMARKER_OUTLIER.value
        measurements_df[Columns.FLAG.value] += np.where((measurements_df[Columns.NUMBER_OF_USABLE_BIOMARKERS.value] <= 1),
                                                 flag_biomarker_outlier, 0)

        biomarker_stat = measurements_df[Columns.NUMBER_OF_USABLE_BIOMARKERS.value].value_counts().sort_index(ascending=True)
        stat_output = ""
        for num_biomarker, num_samples in zip(biomarker_stat.index.values, biomarker_stat.values):
            stat_output += "\t{} usable biomarkers found in {} samples\n".format(num_biomarker, num_samples)
        self.logger.log.info("[Select Biomarkers] - [Sample location: '{}'] - \n"
                             "{}".format(sample_location,  stat_output))

    def report_last_biomarkers_invalid(self, sample_location, measurements_df: pd.DataFrame):
        last_two_measurements = measurements_df.tail(self.report_number_of_biomarker_outlier)
        if last_two_measurements.shape[0] == self.report_number_of_biomarker_outlier:
            for biomarker1, biomarker2 in itertools.combinations(Columns.get_biomarker_columns(), 2):
                valid_ratios = SewageFlag.is_not_flag_set_for_series(last_two_measurements[Columns.get_biomaker_ratio_flag(biomarker1, biomarker2)],
                               SewageFlag.BIOMARKER_RATIO_OUTLIER)
                num_valid_ratios_last_measurements = valid_ratios.sum()
                if num_valid_ratios_last_measurements == 0:  # last samples where marked as outlier
                    biomarker_last_outliers = last_two_measurements[
                        [Columns.DATE.value, biomarker1, biomarker2, biomarker1 + "/" + biomarker2]]
                    self.logger.log.warn("[Report last two biomarker outliers] - [Sample location: '{}'] - "
                                         "The last two measurements for biomarker ratio '{}/{}' are marked as outliers.\n"
                                         "{}".format(sample_location, biomarker1, biomarker2, biomarker_last_outliers))
                    # Todo: handle reporting
        else:
            self.logger.log.debug("[Report last two biomarker outliers] - [Sample location: '{}'] - "
                                  "Less than {} last values available. Skipping... ".format(sample_location, self.report_number_of_biomarker_outlier))