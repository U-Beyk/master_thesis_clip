'''
Assembles the data in the fasta files and the csv file with the predictions into one datframe.

author: U.B.
'''

import pandas as pd
from Bio import SeqIO
from typing import Iterator

class DataframeAssembler:
    '''Class representing the collection of all predictions.'''
    def __init__(self, fasta_filepath: str, predictions_filepath: str) -> None:
        '''Initializes the datafarme object via specified fasta filepath and prediction csv filepath'''
        self.sequences_filepath = fasta_filepath
        self.predicitons_filepath = predictions_filepath
        self.rna_dataframe = self._merge_sequences_predictions_df()

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
    
    def __iter__(self) -> Iterator:
        '''Allows iteration through the dataframe directly in the RNADataFrameAssembler object'''
        return self.rna_dataframe.itertuples(index=False)