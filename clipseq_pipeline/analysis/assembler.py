import pandas as pd
from Bio import SeqIO
from typing import Callable, Self

from .record import CLIPmotiFoldRecord, CLIPmotiCesRecord
from .prediction import CLIPmotiFoldPrediction, CLIPmotiCesPrediction

class CLIPmotiFoldAssembler:
    def __init__(self, fasta_filepath: str, predictions_filepath: str) -> None:
        '''Initializes the datafarme object via specified fasta filepath and prediction csv filepath'''
        self.sequences_filepath = fasta_filepath
        self.predicitons_filepath = predictions_filepath
        self.rna_dataframe = self._merge_sequences_predictions_df()
        self.rna_sequences = self._assemble_rna_sequences()

    def _create_fasta_df(self) -> pd.DataFrame:
        '''Creates a dataframe from the fasta file.'''
        records = []
        for record in SeqIO.parse(self.sequences_filepath, "fasta"):
            records.append({"ID": record.id, "sequence": str(record.seq)})
        return pd.DataFrame(records)
    
    def _create_predictions_df(self) -> pd.DataFrame:
        '''Creates a dataframe from the csv file with the predictions.'''
        df = pd.read_csv(self.predicitons_filepath, sep="\t")
        df = df[df["mfe"] <= 0].reset_index(drop=True)
        return df
    
    def _merge_sequences_predictions_df(self) -> pd.DataFrame:
        '''Merges the dataframe with te sequences and the dataframe with the predictions.'''
        predictions_df = self._create_predictions_df()
        sequences_df = self._create_fasta_df()
        return pd.merge(predictions_df, sequences_df, on="ID")
    
    def _assemble_rna_sequences(self) -> dict[str, CLIPmotiFoldRecord]:
        '''Assembles the RNA sequences with RNAmotiFold records.'''
        rna_sequences: dict[str, CLIPmotiFoldRecord] = {}
        for seq_header, seq_class, seq_mfe, seq_motbracket, sequence in self.rna_dataframe.itertuples(index=False, name=None):
            rna_sequence = rna_sequences.setdefault(seq_header, CLIPmotiFoldRecord.from_header(seq_header, sequence))
            prediction_class = "" if pd.isna(seq_class) else seq_class
            free_energy = seq_mfe / 100
            rna_sequence.add_prediction(CLIPmotiFoldPrediction(prediction_class, free_energy, seq_motbracket, rna_sequence.rel_clip_start_end))
        return rna_sequences
    
    def filter_records(self, predicate: Callable[[CLIPmotiFoldRecord], bool]) -> Self:
        """
        Returns a new assembler containing only records for which `predicate(record)` is True.
        
        This method is functional and chainable — you can call it repeatedly to refine filters.

        Example:
            assembler.filter_records(lambda r: r.mfe_value < -5) \
                     .filter_records(lambda r: "UGCAU" in r.motifs_set)
        """
        new_sequences = {
            key: record
            for key, record in self.rna_sequences.items()
            if predicate(record)
        }
        return self._set_sequences(new_sequences)
    
    def filter_lowest_mfe_predictions(self) -> Self:
        new_sequences = {
            key: record.filter_predictions(0)
            for key, record in self.rna_sequences.items()
        }
        return self._set_sequences(new_sequences)

    def _set_sequences(self, sequences: dict[str, CLIPmotiFoldRecord]) -> Self:
        '''Creates a new instance of the corresponding class, with a specified set of RNA sequences.'''
        new_instance = self.__class__(self.sequences_filepath, self.predicitons_filepath)
        new_instance.rna_sequences = sequences
        return new_instance

# TODO: Refactor class and methods.
class CLIPmotiCesAssembler(CLIPmotiFoldAssembler):
    def filter_records(self, condition: Callable[[CLIPmotiCesRecord], bool] | None = None, mfe_threshold: float | None = None) -> Self:
        """Filter records by condition and/or provided MFE threshold."""
        new_sequences = {}
        for key, record in self.rna_sequences.items():
            filtered_record = record.filter_predictions(mfe_threshold) if mfe_threshold else record
            if condition is None or condition(filtered_record):
                new_sequences[key] = filtered_record
        return self._set_sequences(new_sequences)

    def _set_sequences(self, sequences: dict[str, CLIPmotiCesRecord]) -> Self:
        '''Creates a new instance of the corresponding class, with a specified set of RNA sequences.'''
        new_instance = self.__class__(self.sequences_filepath, self.predicitons_filepath)
        new_instance.rna_sequences = sequences
        return new_instance
    
    def _assemble_rna_sequences(self) -> dict[str, CLIPmotiCesRecord]:
        '''Assembles the RNA sequences with RNAmotiFold records.'''
        rna_sequences: dict[str, CLIPmotiCesRecord] = {}
        for row in self.rna_dataframe.itertuples(index=False):
            rna_sequence = rna_sequences.setdefault(row.ID, CLIPmotiCesRecord.from_header(row.ID, row.sequence))
            prediction_class = "" if pd.isna(row.Class) else row.Class
            free_energy = row.mfe / 100
            rel_start_end = (
                rna_sequence.clip_start - rna_sequence.sequence_start,
                rna_sequence.clip_end - rna_sequence.sequence_start
                )
            rna_sequence.add_prediction(CLIPmotiCesPrediction(prediction_class, free_energy, row.motBracket, rel_start_end))
        return rna_sequences