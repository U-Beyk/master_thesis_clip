import pandas as pd
from collections import Counter

from .record import CLIPmotiFoldRecord, CLIPmotiCesRecord

class CLIPmotiFoldStatistics:

    def __init__(self, rna_sequences: dict[str, CLIPmotiFoldRecord]):
        self.rna_sequences = rna_sequences

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
    def motifs_count(self) -> pd.Series:
        '''Counts all unique motifs for each sequence. Also counts all sequences without any motifs.'''
        motifs_counter = Counter()
        for rna in self.rna_sequences.values():
            motifs_counter.update(rna.motifs_set)
        return pd.Series(motifs_counter).sort_values(ascending=False)
    
    @property
    def clip_overlap_motifs_set(self) -> set[str]:
        motifs = set()
        for rna in self.rna_sequences.values():
            motifs.update(rna.clip_overlap_motifs_set)
        return motifs
    
    @property
    def clip_overlap_motifs_count(self) -> pd.Series:
        motifs_counter = Counter()
        for rna in self.rna_sequences.values():
            motifs_counter.update(rna.clip_overlap_motifs_set)
        return pd.Series(motifs_counter).sort_values(ascending=False)
    
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
    def mfe_motif_distances(self) -> pd.DataFrame:
        '''Creates a dtaframe with the distance to the lowest mfe for each motif.'''
        distance_motifs = []
        for rna_seq in self.rna_sequences.values():
            for motif_data in rna_seq.mfe_motif_distances:
                for motif, value in motif_data.items():
                    row = {"Motifs": motif, "Distance to mfe": value}
                    distance_motifs.append(row)
        return pd.DataFrame(distance_motifs)
    
    @property
    def median_mfe_motif_distances(self) -> pd.Series:
        '''Calculates a series with the median of the distance to the lowest mfe for each motif.'''
        return self.mfe_motif_distances.groupby("Motifs")["Distance to mfe"].median()
    
    @property
    def mfe_clip_motif_distances(self) -> pd.DataFrame:
        distance_motifs = []
        for rna_seq in self.rna_sequences.values():
            for motif_data in rna_seq.mfe_clip_motif_distances:
                for motif, value in motif_data.items():
                    row = {"Motifs": motif, "Distance to mfe": value}
                    distance_motifs.append(row)
        return pd.DataFrame(distance_motifs)

    @property
    def median_mfe_clip_motif_distances(self) -> pd.Series:
        return self.mfe_clip_motif_distances.groupby("Motifs")["Distance to mfe"].median()

# TODO: Change or implement new method, that gets only the potential motifs, that overlap with CLIP data.
class CLIPmotiCesStatistics(CLIPmotiFoldStatistics):

    def __init__(self, rna_sequences: dict[str, CLIPmotiCesRecord]):
        self.rna_sequences = rna_sequences

    @property
    def potential_motifs_count(self) -> dict:
        '''Gets and returns all potential motifs computed via positional abstraction.'''
        motifs_counter = Counter()
        for rna in self.rna_sequences.values():
            motifs_counter.update(rna.potential_motifs_set)
        return dict(motifs_counter.most_common(10))