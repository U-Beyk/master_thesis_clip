from abc import ABC, abstractmethod
from dataclasses import dataclass
from enum import Enum
import json
import os

import pandas as pd
from scipy.stats import chi2_contingency, chisquare

class ReportStrategy(ABC):
    @abstractmethod
    def report(
        self,
        df: pd.DataFrame,
        null_model_df: pd.DataFrame | None = None
    ) -> dict:
        pass

    # TODO: Refactor code.
    def _pivot_motif_counts_by_rbp(self, df: pd.DataFrame) -> dict[str, pd.DataFrame]:
        mot_freq_rbp = {}
        filtered_df = df.xs("motif_frequency", level="format_name")
        for filter_name in filtered_df.index.get_level_values("filter_name").unique():
            rbp_subset = filtered_df.xs(filter_name, level="filter_name")
            combined = []
            if rbp_subset.empty:
                continue
            rbp_data: pd.DataFrame = rbp_subset["data"]
            for rbp_name, sub_df in rbp_data.items():
                if rbp_name == "all_data":
                    continue
                if sub_df is None or sub_df.empty:
                    continue
                temp = sub_df.copy()
                temp["rbp_name"] = rbp_name
                combined.append(temp)
            if not combined:
                continue
            pivot_df = (
                pd.concat(combined, ignore_index=True)
                .pivot(index="Motifs", columns="rbp_name", values="Count")
                .fillna(0)
            )
            mot_freq_rbp[filter_name] = pivot_df
        return mot_freq_rbp

# TODO: Refactor code.
class Chi2ProteinsMotifs(ReportStrategy):
    # H0: Motif occurrence is independent of protein; H1: Motif occurrence depends on protein
    # if p < 0.05 H0 will be rejected
    def report(self, df: pd.DataFrame) -> dict:
        chi2_results = self._chi2_test(df)
        return {
            "description": "Chi-square test of independence across RBPs",
            "results": chi2_results
        }

    def _chi2_test(self, df: pd.DataFrame) -> dict[str, dict[str, float | int]]:
        motif_freq_rbp = self._pivot_motif_counts_by_rbp(df)
        m = len(motif_freq_rbp)
        results = {}
        for filter_, df in motif_freq_rbp.items():
            matrix = df.values
            chi2, p, dof, _ = chi2_contingency(matrix)
            if chi2 == 0:
                continue
            results[filter_] = {
                "chi2": float(chi2),
                "p_value": float(p),
                "dof": int(dof),
            }
        return results
    
# TODO: Add Bonferroni correction AND test it.
class Chi2ProteinsNullComparison(ReportStrategy):
    # p<0.05 distributions are different
    '''
    Plan:
    über filter loopen in beiden dfs -> über proteine loopen in beiden dfs
    -> für jedes Protein chisquare goodness-of-fit test machen

    Dfs richtig formatieren.
    '''
    def report(self, df: pd.DataFrame, df_null: pd.DataFrame) -> str:
        chi2_results = self._chi2_test(df, df_null)
        return {
            "description": "Chi-square test comparing RBP motif distributions to null distributions",
            "results": chi2_results
        }

    def _chi2_test(self, df: pd.DataFrame, df_null: pd.DataFrame) -> dict[str, dict[str, dict[str, float]]]:
        # NOTE: Dictionary looks like this: {filter: {protein: p-value}}
        results = dict()
        filter_dfs = self._pivot_motif_counts_by_rbp(df)
        filter_null_dfs = self._pivot_motif_counts_by_rbp(df_null)

        for filter_, filter_df in filter_dfs.items():
            rbp_null_dfs = filter_null_dfs[filter_]
            if filter_ not in results:
                results[filter_] = dict()
            pvalues = {}
            for rbp in filter_df.columns:
                rbp_df = filter_df[[rbp]]
                rbp_null_df = rbp_null_dfs[[rbp_df.columns[0]]]
                rbp_df, rbp_null_df = self._align_on_motifs(rbp_df, rbp_null_df)
                chi2, pvalue = self._chi_square_protein(rbp_df, rbp_null_df)
                pvalues[rbp] = {"chi2": chi2, "p-value": pvalue}
            results[filter_] = self._apply_bonferroni(pvalues)

        return results
    
    def _pivot_motif_counts_by_rbp(self, df: pd.DataFrame) -> dict[str, pd.DataFrame]:
        mot_freq_rbp = {}
        filtered_df = df.xs("motif_frequency", level="format_name")
        for filter_name in filtered_df.index.get_level_values("filter_name").unique():
            rbp_subset = filtered_df.xs(filter_name, level="filter_name")
            combined = []
            if rbp_subset.empty:
                continue
            rbp_data: pd.DataFrame = rbp_subset["data"]
            for rbp_name, sub_df in rbp_data.items():
                if sub_df is None or sub_df.empty:
                    continue
                temp = sub_df.copy()
                temp["rbp_name"] = rbp_name
                combined.append(temp)
            if not combined:
                continue
            pivot_df = (
                pd.concat(combined, ignore_index=True)
                .pivot(index="Motifs", columns="rbp_name", values="Count")
                .fillna(0)
            )
            mot_freq_rbp[filter_name] = pivot_df
        return mot_freq_rbp

    def _extract_protein_df(self, df: pd.DataFrame) -> dict[str, pd.DataFrame]:
        rbp_dfs = {}
        for rbp in df.columns:
            rbp_dfs[rbp] = df[[rbp]]
        return rbp_dfs
    
    def _align_on_motifs(
            self,
            df: pd.DataFrame,
            df_null: pd.DataFrame
        ) -> tuple[pd.DataFrame, pd.DataFrame]:
        all_motifs = df.index.union(df_null.index)
        df = df.reindex(all_motifs).fillna(0)
        df_null = df_null.reindex(all_motifs).fillna(0)
        return df, df_null
    
    def _chi_square_protein(self, df: pd.DataFrame, null_df: pd.DataFrame) -> tuple[float, float]:
        mask = null_df.iloc[:, 0] > 0
        observed = df.iloc[:, 0][mask]
        expected = null_df.iloc[:, 0][mask]
        expected = expected * (observed.sum() / expected.sum())
        chi, p = chisquare(f_obs=observed, f_exp=expected)
        return chi, p

    def _apply_bonferroni(self, pvalue_dict: dict[str, dict[str, float]]) -> dict[str, dict[str, float]]:
        n_tests = len(pvalue_dict)
        corrected = {}
        for rbp, vals in pvalue_dict.items():
            adjusted_p = min(vals["p-value"] * n_tests, 1.0)
            corrected[rbp] = {
                "chi2": vals["chi2"],
                "p-value": vals["p-value"],
                "p-value_bonferroni": adjusted_p
            }
        return corrected


class ReportModus(Enum):
    SINGLE = "single"
    COMPARATIVE = "comparative"

@dataclass
class RnaReporter:
    report_modus: ReportModus
    report_cls: ReportStrategy
    report_name: str


def write_report(content: dict, filepath: str) -> None:
    os.makedirs(os.path.dirname(filepath), exist_ok=True)
    with open(filepath, "w") as file:
        json.dump(content, file, indent=4)

def apply_statistics(
    formatted_df: pd.DataFrame,
    null_model_df: pd.DataFrame,
    reporters: list[RnaReporter],
    filepath: str
) -> None:
    report_file = f"{filepath}/report.json"
    json_output = dict()
    if formatted_df.empty:
        return 
    for report in reporters:
        if report.report_modus == ReportModus.SINGLE:
            report_content = report.report_cls.report(formatted_df)
        elif report.report_modus == ReportModus.COMPARATIVE:
            if null_model_df.empty:
                continue
            report_content = report.report_cls.report(formatted_df, null_model_df)
        json_output[report.report_name] = report_content
    write_report(json_output, report_file)

RNAMOTIFOLD_REPORTS: list[RnaReporter] = [
    RnaReporter(ReportModus.SINGLE, Chi2ProteinsMotifs(), "chi2_independence",),
    RnaReporter(ReportModus.COMPARATIVE, Chi2ProteinsNullComparison(), "chi2_goodness_of_fit")
]

def report_motifold_df(
        formatted_df: pd.DataFrame,
        null_model_df:pd.DataFrame,
        base_filepath: str
    ) -> None:
    return apply_statistics(formatted_df, null_model_df, RNAMOTIFOLD_REPORTS, base_filepath)


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
    
    data_null = {('motif_frequency', 'mfe_below_zero', 'ms1'): {'M': 87, 'G': 133, 'U': 15, 'T': 98, 'D': 16, 'L': 72},
        ('motif_frequency', 'mfe_below_zero', 'nova3'): {'M': 93, 'G': 74, 'U': 1, 'T': 147,},
        ('motif_frequency', 'mfe_below_zero', 'fmrp'): {'M': 76, 'G': 81, 'U': 13, 'T': 13, 'D': 57, 'L': 41},
        ('motif_frequency', 'mrna_mfe_below_zero', 'ms1'): {'M': 19, 'G': 77, 'U': 162, 'T': 11, 'D': 99, 'L': 64},
        ('motif_frequency', 'mrna_mfe_below_zero', 'nova3'): {'M': 86, 'G': 164, 'U': 119, 'T': 76, 'D': 133, 'L': 58},
        ('motif_frequency', 'mfe_lower_decile', 'ms1'): {'M': 165, 'G': 134, 'U': 91, 'T': 16, 'D': 183, 'L': 5},
        ('motif_frequency', 'mfe_lower_decile', 'nova3'): {'M': 84, 'U': 95, 'T': 173, 'D': 124, 'L': 98},
        ('motif_distances', 'mfe_below_zero', 'ms1'): {'M': 3.77, 'G': 2.10, 'U': 4.64, 'T': 0.84, 'D': 5.32, 'L': 1.17},
        ('motif_distances', 'mfe_below_zero', 'nova3'): {'M': 2.74, 'G': 5.25, 'U': 1.62, 'T': 0.97, 'D': 3.40, 'L': 4.75},
        ('motif_distances', 'mfe_below_zero', 'fmrp'): {'M': 3.56, 'G': 0.86, 'U': 2.47, 'T': 5.46, 'D': 4.08, 'L': 2.67},
        ('motif_distances', 'mrna_mfe_below_zero', 'ms1'): {'G': 3.44, 'U': 0.97, 'T': 4.21, 'D': 3.55, 'L': 2.15},
        ('motif_distances', 'mrna_mfe_below_zero', 'nova3'): {'M': 2.99, 'G': 1.44, 'L': 3.01},
        ('motif_distances', 'mfe_lower_decile', 'ms1'): {'M': 3.11, 'G': 0.99, 'T': 4.88, 'D': 1.20, 'L': 5.02},
        ('motif_distances', 'mfe_lower_decile', 'nova3'): {'M': 4.5, 'U': 0.2, 'T': 3.7, 'D': 5.18, 'L': 1.33}}

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
    df_null = generate_test_df(data_null)
    #TEST BELOW:
    print(Chi2ProteinsNullComparison().report(df, df_null))
    #print(Chi2ProteinsMotifs().report(df))