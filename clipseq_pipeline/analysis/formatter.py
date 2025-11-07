from abc import ABC, abstractmethod
from collections import Counter
from dataclasses import dataclass
from typing import Callable

import pandas as pd

from .transformer import FilteredRbpDf

type FormatFn = Callable[[pd.DataFrame], pd.DataFrame]

@dataclass(frozen=True)
class FormattedRbpDf:
    filter_name: str
    rbp_name: str
    format_name: str
    df: pd.DataFrame

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
    

def apply_formatters(
        filtered_dfs: list[FilteredRbpDf],
        formatters: list[RnaFormatter]
    ) -> list[FormattedRbpDf]:
    formatted_dfs: list[FormattedRbpDf] = []

    for filtered in filtered_dfs:
        for formatter in formatters:
            df = filtered.df
            for fn in formatter.pipeline:
                df = fn(df)
            if df.empty:
                continue
            formatted_dfs.append(
                FormattedRbpDf(
                    filter_name=filtered.filter_name,
                    rbp_name=filtered.rbp_name,
                    format_name=formatter.name,
                    df=df
                )
            )
    return formatted_dfs


RNAMOTIFOLD_FORMATS: list[RnaFormatter] = [
    RnaFormatter("motif_frequency", [MotifFrequencyFormatter().format]),
    RnaFormatter("motif_distances", [MotifDistanceFormatter().format]),
]

def format_motifold_dfs(filtered_dfs: list[FilteredRbpDf]) -> list[FormattedRbpDf]:
    return apply_formatters(filtered_dfs, RNAMOTIFOLD_FORMATS)


RNAMOTICES_FORMATS: list[RnaFormatter] = [
    RnaFormatter("potential_motifs", [PotentialMotifFormatter().format])
]

def format_motices_dfs(filtered_dfs: list[FilteredRbpDf]) -> list[FormattedRbpDf]:
    return apply_formatters(filtered_dfs, RNAMOTICES_FORMATS)