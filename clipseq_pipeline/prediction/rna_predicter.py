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
    filename_types: dict[str, str]
        Dictionary of the different files to use the FASTAs of
        and run the predictions.
    '''

    def __init__(self):
        '''
        Instantiates an RNAPredicter object.
        '''
        self.organisms: dict[str, str] = ORGANISMS_TO_EXAMINE
        self.filename_types = {"rbp_sites": "prediction",
                          "shuffled_rbp_sites": "shuffled_prediction",
                        }
    
    @staticmethod
    def _predict_organism(organism: str, fasta_suffix: str, prediction_suffix: str) -> None:
        '''
        Runs the algorithm to predict the RNA structures of the FASTA file
        of the specified organism.

        Parameters
        ----------
        organism: str
            The string of the organism to predict the structures of.
        fasta_suffix: str
            Suffix for the FASTA file.
        prediciton_suffix: str
            Suffix for the CSV file with the predictions.
        '''
        script_folder = "./RNAmotiFold"
        cmd = [
            "python3",
            "RNAmotiFold.py",
            "-i", f"../data/fasta_files/{organism}/{organism}_{fasta_suffix}.fasta",
            "-o", f"../data/rna_predictions/{organism}/{organism}_{prediction_suffix}.csv",
            "-k", f"{RNA_PREDICTION_K_VALUE}",
            "-w", f"{RNA_PREDICTION_WORKER}",
            "--no_update",
        ]
        subprocess.run(cmd, cwd=script_folder, check=True)

    @staticmethod
    def _delete_old_prediction_file(organism: str, prediciton_suffix: str) -> None:
        '''
        Deletes the old prediction file, so that the predictions won't be appended
        when rerunning the predictions.
        
        Parameters
        ----------
        organism: str
            The string of the organism to delete the prediction file of.
        prediciton_suffix: str
            Suffix for the CSV file with the predictions.
        '''
        file_path = f"./data/rna_predictions/{organism}/{organism}_{prediciton_suffix}.csv"
        if os.path.exists(file_path):
            os.remove(file_path)

    @staticmethod
    def _create_prediction_dir(organism: str) -> None:
        """
        Creates the directory for the predictions of a specified organism.

        Parameters
        ----------
        organism: str
            The string of the organism to create a directory for.
        """
        output_dir = f"./data/rna_predictions/{organism}"
        os.makedirs(output_dir, exist_ok=True)
        
    def run_predictions(self) -> None:
        '''
        Runs the predictions on all specified organisms.
        '''
        for name, organism in self.organisms.items():
            for fasta_name, prediction_name in self.filename_types.items():
                self._create_prediction_dir(organism)
                self._delete_old_prediction_file(organism, prediction_name)
                self._predict_organism(organism, fasta_name, prediction_name)
            print(f"Finished predictions of: {name}!")