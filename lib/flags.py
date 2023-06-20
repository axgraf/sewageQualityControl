# Created by alex at 15.06.23

from enum import Enum

class Flag(Enum):

    COMMENT_NOT_EMPTY = "Comment field not empty"
    BIOMARKER_BELOW_THRESHOLD = "Biomarker values are all empty or below threshold"
    BIOMARKER_DIFFERENCE_TOO_HIGH = "Biomarker ratios are too high compared to last measurements"
    NOT_ENOUGH_BIOMARKERS = "Not enough biomarkers"
