'''
rna_predicter.py
================

Thid module coontains the class running the predictions.

author: U.B.
'''

import os
import subprocess

from ..constants import ORGANISMS_TO_EXAMINE, RNA_PREDICTION_WORKER, RNA_PREDICTION_K_VALUE

class RNAPredicter:
    '''
    Class running the predictions of all specified organisms.

    Attributes
    ----------
    organisms: dict[str, str]
        Dictionary of all organisms to run the predictions of.
    '''

    def __init__(self):
        '''
        Instantiates an RNAPredicter object.
        '''
        self.organisms: dict[str, str] = ORGANISMS_TO_EXAMINE
    
    @staticmethod
    def _predict_organism(organism: str) -> None:
        '''
        Runs the algorithm to predict the RNA structures of the FASTA file
        of the specified organism.

        Parameters
        ----------
        organism: str
            The string of the organism to predict the structures of.
        '''
        script_folder = "./RNAmotiFold"
        cmd = [
            "python3",
            "RNAmotiFold.py",
            "-i", f"../data/fasta_files/{organism}_rbp_sites.fasta",
            "-o", f"../data/rna_predictions/{organism}_prediction.csv",
            "-k", f"{RNA_PREDICTION_K_VALUE}",
            "-w", f"{RNA_PREDICTION_WORKER}",
            "--no_update",
        ]
        subprocess.run(cmd, cwd=script_folder, check=True)

    @staticmethod
    def _delete_old_prediction_file(organism: str) -> None:
        '''
        Deletes the old prediction file, so that the predictions won't be appended
        when rerunning the predictions.
        
        Parameters
        ----------
        organism: str
            The string of the organism to delete the prediction file of.
        '''
        file_path = f"./data/rna_predictions/{organism}_prediction.csv"
        if os.path.exists(file_path):
            os.remove(file_path)
        
    def run_predictions(self) -> None:
        '''
        Runs the predictions on all organisms specified in the config file.
        '''
        for name, organism in self.organisms.items():
            self._delete_old_prediction_file(organism)
            self._predict_organism(organism)
            print(f"Finished predictions of: {name}!")