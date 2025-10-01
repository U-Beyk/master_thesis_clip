'''
rna_predicter.py
~~~~~~~~~~~~~~~~

Classes
-------
RNAPredicter
    Class predicting the RNA-sequences in the FASTA files of all organisms specified in the config file.

author: U.B.
'''

import yaml
import subprocess
import os

class RNAPredicter:
    '''
    Class running the predictions of all organisms specified in the config.yaml.

    Attributes
    ----------
    organisms: dict[str, str]
        Dictionary of all organisms run the RNAmotiFold predictions of.
    '''

    def __init__(self, config_file: str):
        '''
        Instantiates a RNAPredicter object.

        Parameters
        ----------
        config_file: str
            The path to the config file.
        '''
        self.organisms: dict[str, str] = self._extract_organisms(config_file)
    
    @staticmethod
    def _extract_organisms(config_file: str) -> dict[str, str]:
        '''
        Extracts the dictionary containing the organisms to predict the RNA structures of.

        Parameters
        ----------
        config_file: str
            Path to the config file.

        Returns
        -------
        dict[str, str]
            Dictionary of all organisms.
        '''
        with open(config_file, "r") as file:
            content = yaml.safe_load(file)
        return content["organisms"]
    
    @staticmethod
    def _predict_organism(organism: str) -> None:
        '''
        Runs the RNAmotiFold algorithm to predict the RNA structures of the FASTA file
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
            "-k", "10",
            "-w", "12",
            "--no_update",
        ]
        subprocess.run(cmd, cwd=script_folder, check=True)

    @staticmethod
    def _delete_old_prediction_file(organism: str) -> None:
        '''
        Deletes the old prediction file, so that the predictions won't be appended when rerunning pipeline.
        
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
        print("Running predictions!")
        for organism in self.organisms.values():
            self._delete_old_prediction_file(organism)
            self._predict_organism(organism)
        print("Finished predictions!")