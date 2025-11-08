from dataclasses import dataclass

import pandas as pd

@dataclass
class FilterData:
    filter_name: str
    proteins: dict[str, pd.DataFrame]

@dataclass
class AlgorithmData:
    alg_name: str
    filters: dict{str, FilterData}