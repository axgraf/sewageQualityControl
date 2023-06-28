# Created by alex at 15.06.23

from enum import Enum

class Flag(Enum):

    COMMENT_NOT_EMPTY = "Comment field not empty"
    BIOMARKER_BELOW_THRESHOLD = "Biomarker values are all empty or below threshold"
    BIOMARKER_NOT_IN_CONFIDENCE_INTERVAL_99 = "Biomarker ratios outside of 99% confidence interval"
    NOT_ENOUGH_BIOMARKERS_FOR_OUTLIER_DETECTION = "Not enough previous biomarker values for outlier detection"
    BIOMARKER_OUTLIER_FLAG = "At least one biomarker was detected as outlier"
    MIN_BIOMARKER_NUMBER_NOT_REACHED = "Minimal number of biomarkers is not sufficient"


