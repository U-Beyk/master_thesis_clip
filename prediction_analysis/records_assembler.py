'''
Contains all classes that assemble different types of RNA records into a dictionary.

author: U.B.
'''

import pandas as pd
from abc import ABC, abstractmethod
from collections import Counter
from typing import Callable, Self

from prediction_analysis.prediction import CLIPmotiFoldPrediction, CLIPmotiCesPrediction
from prediction_analysis.record import CLIPmotiFoldRecord, CLIPmotiCesRecord
from prediction_analysis.dataframe_assembler import DataframeAssembler
from prediction_analysis.visualization import BarChart, Histogram, ViolinPlot

#TO-DO: Implement the visualization classes taking no dataframe types and instead a dict or count of motifs/motices for consistency in types.
class CLIPRecordsAssembler(ABC):
    def __init__(self, rna_dataframe: DataframeAssembler):
        '''Initializes the assembler via specified dataframe.'''
        self.rna_dataframe = rna_dataframe
        self.rna_sequences = self._assemble_rna_sequences()

    @property
    def unique_rbps(self) -> set[str]:
        unique_rbps = set()
        for rna in self.rna_sequences.values():
            unique_rbps.add(rna.rbp_name)
        return unique_rbps

    @property
    def lowest_mfe_value(self) -> float:
        '''Returns the lowest mfe value as a float.'''
        value = 0
        for record in self.rna_sequences.values():
            if record.mfe_value < value:
                value = record.mfe_value
        return value
    
    @property
    def motifs_set(self) -> set[str]:
        '''Returns a set of unique motifs of all RNA sequences.'''
        motifs = set()
        for rna in self.rna_sequences.values():
            motifs.update(rna.motifs_set)
        return motifs
    
    @property
    def motifs_count(self) -> dict[str, int]:
        '''Counts all unique motifs for each sequence. Also counts all sequences without any motifs.'''
        motifs_counter = Counter()
        for rna in self.rna_sequences.values():
            motifs_counter.update(rna.motifs_set)
        return dict(motifs_counter)
    
    @property
    def clip_overlap_motifs_set(self) -> set[str]:
        motifs = set()
        for rna in self.rna_sequences.values():
            motifs.update(rna.clip_overlap_motifs_set)
        return motifs
    
    @property
    def clip_overlap_motifs_count(self) -> dict[str, int]:
        motifs_counter = Counter()
        for rna in self.rna_sequences.values():
            motifs_counter.update(rna.clip_overlap_motifs_set)
        return dict(motifs_counter)
    
    @property
    def sequence_number(self) -> int:
        '''Gets and returns the number of seuquences.'''
        return len(self.rna_sequences)
    
    @property
    def prediction_number(self) -> int:
        '''Counts and returns the total number of predictions.'''
        number = 0
        for sequence in self.rna_sequences.values():
            number += sequence.prediction_number
        return number
    
    @property
    def distance_to_lowest_all_motifs(self) -> pd.DataFrame:
        '''Creates a dtaframe with the distance to the lowest mfe for each motif.'''
        distance_motifs = []
        for rna_seq in self.rna_sequences.values():
            for motif_data in rna_seq.mfe_motif_distances:
                for motif, value in motif_data.items():
                    row = {"Motifs": motif, "Distance to mfe": value}
                    distance_motifs.append(row)
        return pd.DataFrame(distance_motifs)
    
    @property
    def median_distance_to_lowest_all_motifs(self) -> pd.Series:
        '''Calculates a series with the median of the distance to the lowest mfe for each motif.'''
        return self.distance_to_lowest_all_motifs.groupby("Motifs")["Distance to mfe"].median()
    
    @property
    def distance_to_lowest_all_clip_motifs(self) -> pd.DataFrame:
        distance_motifs = []
        for rna_seq in self.rna_sequences.values():
            for motif_data in rna_seq.mfe_clip_motif_distances:
                for motif, value in motif_data.items():
                    row = {"Motifs": motif, "Distance to mfe": value}
                    distance_motifs.append(row)
        return pd.DataFrame(distance_motifs)

    @property
    def median_distance_to_lowest_all_clip_motifs(self) -> pd.Series:
        return self.distance_to_lowest_all_clip_motifs.groupby("Motifs")["Distance to mfe"].median()
    
    def filter_records(self, condition: Callable[[CLIPmotiFoldRecord], bool] | None = None, mfe_range: bool = False) -> Self:
        '''Filters the records by a specified condition and optionally gets the predictions in a specified mfe range.
        Returns a new instance of the corresponding class.'''
        mfe_threshold = self.lowest_mfe_value * 0.1 if mfe_range else None
        new_sequences = {}
        for key, record in self.rna_sequences.items():
            filtered_record = record.filter_predictions(mfe_threshold) if mfe_threshold else record
            if condition is None or condition(filtered_record):
                new_sequences[key] = filtered_record
        return self._set_sequences(new_sequences)

    def _set_sequences(self, sequences: dict[str, CLIPmotiFoldRecord]) -> Self:
        '''Creates a new instance of the corresponding class, with a specified set of RNA sequences.'''
        new_instance = self.__class__(self.rna_dataframe)
        new_instance.rna_sequences = sequences
        return new_instance
    
    @abstractmethod
    def _assemble_rna_sequences(self) -> dict[str, CLIPmotiFoldRecord]:
        pass

    def __repr__(self) -> str:
        '''Represents the object and its details as a string.'''
        lines = [f"{key}: {value}" for key, value in self.rna_sequences.items()]
        return "\n".join(lines)


class CLIPmotiFoldRecordsAssembler(CLIPRecordsAssembler):
    '''Generic class for an RNA record assembler.'''

    # TO-DO: Make IDE recognize the datatypes.
    def _assemble_rna_sequences(self) -> dict[str, CLIPmotiFoldRecord]:
        '''Assembles the RNA sequences with RNAmotiFold records.'''
        rna_sequences: dict[str, CLIPmotiFoldRecord] = {}
        for row in self.rna_dataframe:
            rna_sequence = rna_sequences.setdefault(row.ID, CLIPmotiFoldRecord.from_header(row.ID, row.sequence))
            prediction_class = "" if pd.isna(row.Class) else row.Class
            free_energy = row.mfe / 100
            rel_start_end = (
                rna_sequence.clip_start - rna_sequence.sequence_start,
                rna_sequence.clip_end - rna_sequence.sequence_start
                )
            rna_sequence.add_prediction(CLIPmotiFoldPrediction(prediction_class, free_energy, row.motBracket, rel_start_end))
        return rna_sequences
    
    def visualize_as_barchart(self, description: str, filepath: str) -> None:
        '''Visualizes the data as a barchart.'''
        barchart = BarChart(
            description,
            self.clip_overlap_motifs_count,
            self.sequence_number
        )
        barchart.save_plot(filepath)

    def visualize_as_violinplot(self, description: str, filepath: str) -> None:
        '''Visualizes the data as a violinplot.'''
        plot_data = self.distance_to_lowest_all_clip_motifs
        violinplot = ViolinPlot(
            description,
            plot_data,
            self.prediction_number
        )
        violinplot.save_plot(filepath)
        pass

class CLIPmotiCesRecordsAssembler(CLIPRecordsAssembler):
    '''Generic class for RNAHeliCes/RNAmotiCes record assemblies.'''

    # TO-DO: Adjust or implement potential motif count to CLIP-motif overlap.
    @property
    def potential_motifs_count(self) -> dict:
        '''Gets and returns all potential motifs computed via positional abstraction.'''
        motifs_counter = Counter()
        for rna in self.rna_sequences.values():
            motifs_counter.update(rna.potential_motifs_set)
        return dict(motifs_counter.most_common(10))
    
    # TO-DO: Make IDE recognize the datatypes.
    def _assemble_rna_sequences(self) -> dict[str, CLIPmotiCesRecord]:
        '''Assembles the RNA sequences with RNAHeliCes/RNAmotiCes records.'''
        rna_sequences: dict[str, CLIPmotiCesRecord] = {}
        for row in self.rna_dataframe:
            rna_sequence = rna_sequences.setdefault(row.ID, CLIPmotiCesRecord.from_header(row.ID, row.sequence))
            prediction_class = "" if pd.isna(row.Class) else row.Class
            free_energy = row.mfe / 100
            rel_start_end = (
                rna_sequence.clip_start - rna_sequence.sequence_start,
                rna_sequence.clip_end - rna_sequence.sequence_start
                )
            rna_sequence.add_prediction(CLIPmotiCesPrediction(prediction_class, free_energy, row.motBracket, rel_start_end))
        return rna_sequences
    
    # TO-DO. Adjust method to CLIP-motif overlap.
    def visualize_as_histogram(self, description: str, filepath: str) -> None:
        '''Visualizes the data as a histogram.'''
        histogram = Histogram(
            description,
            self.potential_motifs_count,
            self.sequence_number
        )
        histogram.save_plot(filepath)