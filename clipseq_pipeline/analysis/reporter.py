from dataclasses import dataclass
from functools import reduce
import pandas as pd

import numpy as np
from scipy.stats import chi2_contingency

from .formatter import FormattedRbpDf

@dataclass
class RnaReporter:
    rbp_name: str
    filter_name: str
    format_name: str
    report_name: str
    df: pd.DataFrame


def generate_motif_rbp_df(dfs: list[FormattedRbpDf]) -> dict[str, pd.DataFrame]:
    # Collect all unique filter_names
    filter_names = set(rbp.filter_name for rbp in dfs)

    result = {}

    for f_name in filter_names:
        # Filter dfs for this filter_name and motif_frequency
        motif_dfs = []
        for rbp in dfs:
            if rbp.format_name == "motif_frequency" and rbp.rbp_name != "ALL_DATA" and rbp.filter_name == f_name:
                df = rbp.df.copy()
                if not {"Motifs", "Count"}.issubset(df.columns):
                    raise ValueError(f"{rbp.rbp_name} DataFrame must have 'Motifs' and 'Count' columns.")
                
                # Rename 'Count' → protein name
                df = df.rename(columns={"Count": rbp.rbp_name})
                motif_dfs.append(df)
        
        if not motif_dfs:
            continue  # skip empty groups

        # Merge all DataFrames for this filter_name
        merged_df = reduce(lambda left, right: pd.merge(left, right, on="Motifs", how="outer"), motif_dfs)
        merged_df = merged_df.fillna(0)
        merged_df.iloc[:, 1:] = merged_df.iloc[:, 1:].astype(int)

        result[f_name] = merged_df

    return result

# H0: Motif occurrence is independent of protein; H1: Motif occurrence depends on protein
# if p < 0.05 H0 will be rejected
def chi2_proteins_motifs(dfs: list[FormattedRbpDf]) -> tuple[float, float, int, np.ndarray]:
    # First, generate merged motif × protein DataFrames per filter
    merged_dict = generate_motif_rbp_df(dfs)

    results = {}

    for f_name, df in merged_dict.items():
        # Convert DataFrame to numeric matrix for chi-square test
        # Rows = motifs, Columns = proteins
        matrix = df.iloc[:, 1:].to_numpy(dtype=int)

        # Run chi-square test of independence
        chi2, p, dof, expected = chi2_contingency(matrix)

        # Store results in dict keyed by filter name
        results[f_name] = {"p-value": p}

    return results

# TODO: Add adhoc test?