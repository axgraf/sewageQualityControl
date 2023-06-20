# Created by alex at 20.06.23
from .sewage import SewageSample
from typing import Dict, List
import pandas as pd


def convert_sample_list2pandas(measurements: List[SewageSample]):
    table = pd.DataFrame.from_records([vars(s) for s in measurements])
    return table
