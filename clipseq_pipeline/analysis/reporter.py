from abc import ABC, abstractmethod
from dataclasses import dataclass
import os
from typing import Callable

import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency

type ReportFn = Callable[[pd.DataFrame], None]

@dataclass
class RnaReporter:
    report_name: str
    report_fn: ReportFn

class ReportStrategy(ABC):
    @abstractmethod
    def report(self, df: pd.DataFrame, filepath: str) -> str:
        pass

# TODO: Implement class.
class Chi2ProteinsMotifs(ReportStrategy):
    # H0: Motif occurrence is independent of protein; H1: Motif occurrence depends on protein
    # if p < 0.05 H0 will be rejected
    def report(self, df: pd.DataFrame) -> str:
        chi2_results = self._chi2_test(df)
        report = "## $\u03C7^2$-Test\n\n$\u03C7^2$-test between the rbps and motif for following filters:\n\n"
        for filter_, result in chi2_results.items():
            filter_report = f"{filter_}: p-value: {result["p-value"]}; chi2: {result["chi2"]}\n\n"
            report += filter_report
        return report

    def _chi2_test(self, df: pd.DataFrame) -> dict[str, dict[str, float | int]]:
        motif_freq_rbp = self._pivot_motif_counts_by_rbp(df)
        results = {}
        for filter_, df in motif_freq_rbp.items():
            matrix = df.values
            chi2, p, dof, _ = chi2_contingency(matrix)
            results[filter_] = {"chi2": chi2, "p-value": p, "dof": dof}
        return results

    def _pivot_motif_counts_by_rbp(self, df: pd.DataFrame) -> dict[str, pd.DataFrame]:
        mot_freq_rbp = {}
        filtered_df = df.xs("motif_frequency", level="format_name")
        for filter_name in filtered_df.index.get_level_values("filter_name").unique():
            rbp_subset = filtered_df.xs(filter_name, level="filter_name")
            combined = []
            rbp_data: pd.DataFrame = rbp_subset["data"]
            for rbp_name, sub_df in rbp_data.items():
                if rbp_name == "all_data":
                    continue
                temp = sub_df.copy()
                temp["rbp_name"] = rbp_name
                combined.append(temp)
            pivot_df = (
                pd.concat(combined, ignore_index=True)
                .pivot(index="Motifs", columns="rbp_name", values="Count")
                .fillna(0)
            )
            mot_freq_rbp[filter_name] = pivot_df
        return mot_freq_rbp

# TODO: Add adhoc test?

def write_report(content: str, filepath: str) -> None:
    os.makedirs(os.path.dirname(filepath), exist_ok=True)
    with open(filepath, "w") as file:
        file.write(content)

# TODO: Fix method.
def apply_statistics(
    formatted_df: pd.DataFrame,
    reporters: list[ReportStrategy],
    filepath: str
) -> None:
    if formatted_df.empty:
        return 
    
    for report_fn in reporters:
        report_content = report_fn(formatted_df)
        write_report(report_content, f"{filepath}/report.md")

RNAMOTIFOLD_REPORTS: list[ReportStrategy] = [
    Chi2ProteinsMotifs().report,
]

def report_motifold_df(formatted_df: pd.DataFrame, base_filepath: str) -> None:
    return apply_statistics(formatted_df, RNAMOTIFOLD_REPORTS, base_filepath)


if __name__ == "__main__":
    data = {('motif_frequency', 'mfe_below_zero', 'ms1'): {'C': 148, 'M': 87, 'G': 133, 'U': 115, 'T': 98, 'D': 166, 'L': 72},
        ('motif_frequency', 'mfe_below_zero', 'nova3'): {'C': 102, 'M': 93, 'G': 74, 'U': 185, 'T': 147,},
        ('motif_frequency', 'mfe_below_zero', 'fmrp'): {'C': 90, 'M': 176, 'G': 81, 'U': 130, 'T': 113, 'D': 57, 'L': 41},
        ('motif_frequency', 'mrna_mfe_below_zero', 'ms1'): {'C': 88, 'M': 139, 'G': 77, 'U': 162, 'T': 121, 'D': 99, 'L': 64},
        ('motif_frequency', 'mrna_mfe_below_zero', 'nova3'): {'M': 86, 'G': 164, 'U': 119, 'T': 76, 'D': 133, 'L': 58},
        ('motif_frequency', 'mfe_lower_decile', 'ms1'): {'C': 78, 'M': 165, 'G': 134, 'U': 91, 'T': 156, 'D': 183, 'L': 55},
        ('motif_frequency', 'mfe_lower_decile', 'nova3'): {'C': 112, 'M': 184, 'U': 95, 'T': 173, 'D': 124, 'L': 98},
        ('motif_distances', 'mfe_below_zero', 'ms1'): {'C': 1.83, 'M': 3.77, 'G': 2.10, 'U': 4.64, 'T': 0.84, 'D': 5.32, 'L': 1.17},
        ('motif_distances', 'mfe_below_zero', 'nova3'): {'M': 2.74, 'G': 5.25, 'U': 1.62, 'T': 0.97, 'D': 3.40, 'L': 4.75},
        ('motif_distances', 'mfe_below_zero', 'fmrp'): {'C': 4.31, 'M': 3.56, 'G': 0.86, 'U': 2.47, 'T': 5.46, 'D': 4.08, 'L': 2.67},
        ('motif_distances', 'mrna_mfe_below_zero', 'ms1'): {'C': 2.10, 'G': 3.44, 'U': 0.97, 'T': 4.21, 'D': 3.55, 'L': 2.15},
        ('motif_distances', 'mrna_mfe_below_zero', 'nova3'): {'C': 0.87, 'M': 2.99, 'G': 1.44, 'L': 3.01},
        ('motif_distances', 'mfe_lower_decile', 'ms1'): {'C': 1.55, 'M': 3.11, 'G': 0.99, 'T': 4.88, 'D': 1.20, 'L': 5.02},
        ('motif_distances', 'mfe_lower_decile', 'nova3'): {'C': 2.41, 'M': 4.05, 'U': 0.92, 'T': 3.67, 'D': 5.18, 'L': 1.33}}

    def generate_test_df(data: dict):
        records = []
        for (format_name, filter_name, rbp_name), motifs in data.items():
            value_col = "Count" if "frequency" in format_name else "Mfe distances [kcal/mol]"
            df_inner = pd.DataFrame({
                "Motifs": list(motifs.keys()),
                value_col: list(motifs.values())
            })
            records.append((format_name, filter_name, rbp_name, df_inner))
        result_df = pd.DataFrame(
            records, columns=["format_name", "filter_name", "rbp_name", "data"]
        ).set_index(["format_name", "filter_name", "rbp_name"])
        return result_df

    df = generate_test_df(data)
    #TEST BELOW:
    report_motifold_df(df, ".")