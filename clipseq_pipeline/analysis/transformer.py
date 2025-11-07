from dataclasses import dataclass
from typing import Callable

import pandas as pd

type FilterFn = Callable[[pd.DataFrame], pd.DataFrame]

@dataclass(frozen=True)
class RnaFilter:
    name: str
    pipeline: list[FilterFn]

@dataclass(frozen=True)
class FilteredRbpDf:
    filter_name: str
    rbp_name: str
    df: pd.DataFrame

def filter_mfe_below_zero(df: pd.DataFrame) -> pd.DataFrame:
    return df[df["mfe"] < 0]

def filter_mfe_minimum(df: pd.DataFrame) -> pd.DataFrame:
    return df[df["mfe"] == df.groupby("sequence_id")["mfe"].transform("min")]

def filter_mfe_lower_decile(df: pd.DataFrame) -> pd.DataFrame:
    mfe_range = abs(df["mfe"].min() * 0.1)
    return df[df["mfe_distance"] < mfe_range]

def filter_by_mrna(df: pd.DataFrame) -> pd.DataFrame:
    return df[df["feature_types"].str.contains("mRNA")]

def filter_by_lncrna(df: pd.DataFrame) -> pd.DataFrame:
    return df[df["feature_types"].str.contains("lnc_RNA")]

def filter_by_rrna(df: pd.DataFrame) -> pd.DataFrame:
    return df[df["feature_types"].str.contains("rRNA")]

def filter_by_trna(df: pd.DataFrame) -> pd.DataFrame:
    return df[df["feature_types"].str.contains("tRNA")]

def filter_by_snrna(df: pd.DataFrame) -> pd.DataFrame:
    return df[df["feature_types"].str.contains("snRNA")]

def filter_by_mirna(df: pd.DataFrame) -> pd.DataFrame:
    return df[df["feature_types"].str.contains("miRNA")]


def build_rbp_dfs(df: pd.DataFrame) -> dict[str, pd.DataFrame]:
    dict_dfs = {"ALL_DATA": df.copy()}
    for rbp in df["rbp_name"].unique():
        rbp_df: pd.DataFrame = df[df["rbp_name"] == rbp]
        dict_dfs[rbp] = rbp_df.copy()
    return dict_dfs

def apply_pipeline(df: pd.DataFrame, pipeline: list[FilterFn]) -> pd.DataFrame:
    for fn in pipeline:
        df = fn(df)
    return df

def apply_rna_filters(df: pd.DataFrame, filters: list[RnaFilter]) -> list[FilteredRbpDf]:
    rbp_dfs = build_rbp_dfs(df)
    filtered_dfs: list[FilteredRbpDf] = []

    for rna_filter in filters:
        for rbp_name, rbp_df in rbp_dfs.items():
            filtered_df = apply_pipeline(rbp_df, rna_filter.pipeline)
            if filtered_df.empty:
                continue
            filtered_dfs.append(FilteredRbpDf(rna_filter.name, rbp_name, filtered_df))

    return filtered_dfs

RNAMOTIFOLD_FILTERS: list[RnaFilter] = [
    RnaFilter("mfe_below0", [filter_mfe_below_zero]),
    RnaFilter("mfe_minimum", [filter_mfe_minimum]),
    RnaFilter("mfe_lower_decile", [filter_mfe_lower_decile]),
    RnaFilter("mrna_mfe_below0", [filter_mfe_below_zero, filter_by_mrna]),
    RnaFilter("lncrna_mfe_below0", [filter_mfe_below_zero, filter_by_lncrna]),
    RnaFilter("rrna_mfe_below0", [filter_mfe_below_zero, filter_by_rrna]),
    RnaFilter("snrna_mfe_below0", [filter_mfe_below_zero, filter_by_snrna]),
    RnaFilter("trna_mfe_below0", [filter_mfe_below_zero, filter_by_trna]),
    RnaFilter("mirna_mfe_below0", [filter_mfe_below_zero, filter_by_mirna])
]

def trans_motifold_dfs(df: pd.DataFrame) -> list[FilteredRbpDf]:
    return apply_rna_filters(df, RNAMOTIFOLD_FILTERS)

RNAMOTICES_FILTERS: list[RnaFilter] = [
    # TODO: Add filters.
]

# TODO: Implement method.
def trans_motices_dfs(df:pd.DataFrame) -> list[FilteredRbpDf]:
    raise NotImplementedError("trans_motices_dfs method is not implemented yet.")