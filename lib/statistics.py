# Created by alex at 14.07.23


class SewageStat:

    def __init__(self):
        self.sample_location = ""
        self.total_samples = 0
        self.stat_dict = dict()
        self.biomarker_dict = dict()
        self.biomarker_ratio_outliers = dict()
        self.surrogate_virus_outliers = dict()
        self.sewage_flow_outliers = dict()
        self.water_quality_outliers = dict()
        self.normalization_outliers = dict()
        self.outliers = dict()

    def __reset(self):
        self.stat_dict = dict()
        self.biomarker_dict = dict()
        self.biomarker_ratio_outliers = dict()
        self.surrogate_virus_outliers = dict()
        self.sewage_flow_outliers = dict()
        self.water_quality_outliers = dict()
        self.normalization_outliers = dict()
        self.outliers = dict()

    def set_sample_location_and_total_number(self, sample_location, total_samples_number):
        self.sample_location = sample_location
        self.total_samples = total_samples_number
        self.__reset()

    def add_comment_not_empty(self):
        msg = 'Comment not empty'
        self.stat_dict.setdefault(msg, 0)
        self.stat_dict[msg] += 1

    def add_mean_sewage_flow_empty(self):
        msg = 'Mean sewage flow empty'
        self.stat_dict.setdefault(msg, 0)
        self.stat_dict[msg] += 1

    def add_biomarker_below_threshold_or_empty(self, biomarker):
        self.biomarker_dict.setdefault(biomarker, 0)
        self.biomarker_dict[biomarker] += 1

    def add_biomarker_ratio_outlier(self, biomarker1, biomarker2, status):
        ratio_name = biomarker1 + "/" + biomarker2
        if not ratio_name in self.biomarker_ratio_outliers:
            self.biomarker_ratio_outliers[ratio_name] = dict()
        self.biomarker_ratio_outliers[ratio_name].setdefault(status, 0)
        self.biomarker_ratio_outliers[ratio_name][status] += 1

    def add_surrogate_virus_outlier(self, sVirus, status):
        if not sVirus in self.surrogate_virus_outliers:
            self.surrogate_virus_outliers[sVirus] = dict()
        self.surrogate_virus_outliers[sVirus].setdefault(status, 0)
        self.surrogate_virus_outliers[sVirus][status] += 1

    def add_sewage_flow_outlier(self, status):
        self.sewage_flow_outliers.setdefault(status, 0)
        self.sewage_flow_outliers[status] += 1

    def add_water_quality_outlier(self, qual_type, status):
        if not qual_type in self.water_quality_outliers:
            self.water_quality_outliers[qual_type] = dict()
        self.water_quality_outliers[qual_type].setdefault(status, 0)
        self.water_quality_outliers[qual_type][status] += 1

    def add_reproduction_factor_outlier(self, status):
        self.normalization_outliers.setdefault(status, 0)
        self.normalization_outliers[status] += 1

    def add_outliers(self, type):
        self.outliers.setdefault(type, 0)
        self.outliers[type] += 1

    def print_statistics(self):
        stats = "{}:\n".format(self.sample_location)
        for key, msg in self.stat_dict.items():
            stats += "{}:\t{}\n".format(key, msg)
        stats += "\nBiomarker values below threshold/empty:\n"
        for biomarker, count in self.biomarker_dict.items():
            stats += "\t{}:\t{}\n".format(biomarker, count)
        stats += "Biomarker ratio status:\n"
        for biomarker_ratio, status_dict in self.biomarker_ratio_outliers.items():
            stats += "\t{}:\n".format(biomarker_ratio)
            for status, count in status_dict.items():
                stats +="\t\t{}:\t{}\n".format(status, count)
        stats += "Sewage flow status:\n"
        for status, count in self.sewage_flow_outliers.items():
            stats += "\t{}:\t{}\n".format(status, count)
        stats += "Water quality status:\n"
        for qual_type, status_dict in self.water_quality_outliers.items():
            stats += "\t{}:\n".format(qual_type)
            for status, count in status_dict.items():
                stats += "\t\t{}:\t{}\n".format(status, count)
        stats += "Outlier reasons:\n"
        for outlier_type, count in self.outliers.items():
            stats += "\t{}:\t{}\n".format(outlier_type, count)
        return stats