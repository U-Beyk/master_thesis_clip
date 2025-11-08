from dataclasses import dataclass
from typing import Callable

import pandas as pd

type FilterFn = Callable[[pd.DataFrame], pd.DataFrame]

@dataclass(frozen=True)
class RnaFilter:
    name: str
    pipeline: list[FilterFn]

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
    dict_dfs = {"all_data": df.copy()}
    for rbp in df["rbp_name"].unique():
        rbp_df: pd.DataFrame = df[df["rbp_name"] == rbp]
        rbp: str = rbp.lower()
        dict_dfs[rbp] = rbp_df.copy()
    return dict_dfs

def apply_pipeline(df: pd.DataFrame, pipeline: list[FilterFn]) -> pd.DataFrame:
    for fn in pipeline:
        df = fn(df)
    return df

# TODO: Refactor code.
def apply_rna_filters_multiindex(df: pd.DataFrame, filters: list[RnaFilter]) -> pd.DataFrame:
    """
    Apply RNA filters to the RBP DataFrames and return a MultiIndex DataFrame.
    
    MultiIndex levels:
    - filter_name
    - rbp_name
    - original row index of filtered df
    """
    rbp_dataframes = build_rbp_dfs(df)
    filtered_results = []

    for rna_filter in filters:
        for rbp_name, rbp_df in rbp_dataframes.items():
            filtered_df = apply_pipeline(rbp_df, rna_filter.pipeline)
            if filtered_df.empty:
                continue

            filtered_df = filtered_df.copy()
            filtered_df["filter_name"] = rna_filter.name
            filtered_df["rbp_name"] = rbp_name

            filtered_results.append(filtered_df)

    if not filtered_results:
        return pd.DataFrame(columns=df.columns).set_index(["filter_name", "rbp_name"])

    combined_df = pd.concat(filtered_results)
    combined_df.set_index(["filter_name", "rbp_name"], inplace=True)
    return combined_df

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

def trans_motifold_dfs(df: pd.DataFrame) -> pd.DataFrame:
    return apply_rna_filters_multiindex(df, RNAMOTIFOLD_FILTERS)

RNAMOTICES_FILTERS: list[RnaFilter] = [
    # TODO: Add filters.
]

# TODO: Implement method.
def trans_motices_dfs(df:pd.DataFrame) -> pd.DataFrame:
    raise NotImplementedError("trans_motices_dfs method is not implemented yet.")