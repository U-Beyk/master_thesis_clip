from abc import ABC, abstractmethod
from collections import Counter
from dataclasses import dataclass
from typing import Callable

import pandas as pd

type FormatFn = Callable[[pd.DataFrame], pd.DataFrame]

@dataclass(frozen=True)
class RnaFormatter:
    name: str
    pipeline: list[FormatFn]

class DataFormatterStrategy(ABC):
    @abstractmethod
    def format(self, df: pd.DataFrame) -> pd.DataFrame:
        """Transform the given DataFrame for a specific plot."""
        pass

class MotifFrequencyFormatter(DataFormatterStrategy):
    def format(self, df: pd.DataFrame) -> pd.DataFrame:
        sequence_motifs = (
            df.groupby("sequence_id")["motifs_clip_window"]
            .apply(lambda motifs: "".join(dict.fromkeys("".join(motifs))))
        )
        barplot_df = self._count_motifs(sequence_motifs)
        return barplot_df
    
    def _count_motifs(self, sequence_motifs: pd.Series) -> pd.DataFrame:
        motif_counts = Counter()
        for motif_str in sequence_motifs:
            motifs_to_count = motif_str or ["No motif"]
            motif_counts.update(motifs_to_count)

        return pd.DataFrame(
            motif_counts.items(), columns=["Motifs", "Count"]
        ).sort_values(by="Count", ascending=False, ignore_index=True)

class MotifDistanceFormatter(DataFormatterStrategy):
    def format(self, df: pd.DataFrame) -> pd.DataFrame:
        motifs_distances = (
            df.groupby("sequence_id")
            .apply(lambda row: list(zip(row["motifs_clip_window"], row["mfe_distance"])))
            .reset_index(name="motif_to_mfe")
            .sort_values("sequence_id")
            .reset_index(drop=True)
        )
        motifs_distances["motif_to_mfe"] = (
            motifs_distances["motif_to_mfe"]
            .apply(self._adjust_motifs)
            .apply(self._extract_min_distances)
        )
        violinplot_df = self._generate_motif_distance_df(motifs_distances)
        return violinplot_df
    
    def _adjust_motifs(self, motifs_list: list[tuple[str, float]]) -> list[tuple[str, float]]:
        adjusted_motifs_list = []
        for motifs, distance in motifs_list:
            distance = round(distance, 1) 
            if len(motifs) > 1:
                adjusted_motifs_list.extend((motif, distance) for motif in motifs)
            elif not motifs:
                adjusted_motifs_list.append(("No motif", distance))
            else:
                adjusted_motifs_list.append((motifs, distance))
        return adjusted_motifs_list
    
    def _extract_min_distances(self, motifs_list: list[tuple[str, float]]) -> list[tuple[str, float]]:
        min_distance_dict = {}
        for motif, distance in motifs_list:
            if motif not in min_distance_dict or distance < min_distance_dict[motif]:
                min_distance_dict[motif] = distance
        
        min_distance_motifs = [(motif, distance) for motif, distance in min_distance_dict.items()]
        return min_distance_motifs
    
    def _generate_motif_distance_df(self, df: pd.DataFrame) -> pd.DataFrame:
        motif_to_distance = []
        for motif_pairs in df["motif_to_mfe"]:
            for motif, distance in motif_pairs:
                motif_to_distance.append((motif, distance))
        return pd.DataFrame(motif_to_distance, columns=["Motifs", "Mfe distances [kcal/mol]"])

# TODO: Implement class.
class PotentialMotifFormatter(DataFormatterStrategy):
    def format(df: pd.DataFrame) -> pd.DataFrame:
        raise NotImplementedError("PotentialMotifFormatter Class and format method are not implemented yet.")
    
# TODO: Refactor code.
def apply_formatters(
        filtered_dfs: pd.DataFrame,
        formatters: list[RnaFormatter]
    ) -> pd.DataFrame:
    """
    Apply a list of RnaFormatter pipelines to each (filter_name, rbp_name) group
    in a MultiIndex DataFrame and return a new DataFrame with a 3-level MultiIndex:
        format_name -> filter_name -> rbp_name
    Each row contains the resulting formatted DataFrame for that group.
    """
    records = []

    grouped = filtered_dfs.groupby(level=["filter_name", "rbp_name"])

    for (filter_name, rbp_name), group_df in grouped:
        for formatter in formatters:
            df = group_df.copy()
            for fn in formatter.pipeline:
                df = fn(df)
            if df.empty:
                continue
            records.append((formatter.name, filter_name, rbp_name, df))

    if not records:
        return pd.DataFrame(
            columns=["format_name", "filter_name", "rbp_name", "data"]
        ).set_index(["format_name", "filter_name", "rbp_name"])

    formatted_df = pd.DataFrame(
        records, columns=["format_name", "filter_name", "rbp_name", "data"]
    ).set_index(["format_name", "filter_name", "rbp_name"])

    return formatted_df


RNAMOTIFOLD_FORMATS: list[RnaFormatter] = [
    RnaFormatter("motif_frequency", [MotifFrequencyFormatter().format]),
    RnaFormatter("motif_distances", [MotifDistanceFormatter().format]),
]

def format_motifold_df(filtered_df: pd.DataFrame) -> pd.DataFrame:
    return apply_formatters(filtered_df, RNAMOTIFOLD_FORMATS)


RNAMOTICES_FORMATS: list[RnaFormatter] = [
    RnaFormatter("potential_motifs", [PotentialMotifFormatter().format])
]

def format_motices_df(filtered_df: pd.DataFrame) -> pd.DataFrame:
    return apply_formatters(filtered_df, RNAMOTICES_FORMATS)