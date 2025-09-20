"""
fasta_builder.py
~~~~~~~~~~~~~~

Contains the class to build FASTA files.

Classes
-------
FastaBuilder
    Class building the FASTA file of organisms.

author: U.B.
"""

import yaml
from concurrent.futures import ProcessPoolExecutor, as_completed

from processing.rbp_site_generator import RbpSiteGenerator

class FastaBuilder:
    '''
    Class creating the FASTA file of specific organisms.

    Attributes
    ----------
    organisms: dict[str, str]
        Dictionary of all organisms to create FASTAs of.
    '''

    def __init__(self, config_file: str):
        '''
        Instantiates a FastaBuilder object.

        Parameters
        ----------
        config_file: str
            Path to the config.yaml with the dictionary of the organisms.
        '''
        self.organisms: dict[str, str] = self._extract_organisms(config_file)
    
    @staticmethod
    def _extract_organisms(config_file: str) -> dict[str, str]:
        '''
        Extracts the dictionary containing the organisms to build FASTAs of.

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
    def _process_organism(organism: str) -> None:
        '''
        Creates a RbpSiteGenerator instance of the specified organism and writes a FASTA file
        of the CLIP data and its sequences.

        Parameters
        ----------
        organism: str
            The string of the organism to process.
        '''
        rbpsite_generator = RbpSiteGenerator(organism)
        rbpsite_generator.write_fasta()

    def build_fasta_organisms(self) -> None:
        '''
        Builds FASTA files of all organisms specified in the config file.
        The FASTA files are build separately each with its own worker.
        '''
        print("Starting to build FASTA files!")
        with ProcessPoolExecutor(max_workers=7) as executor:
            futures = [
                executor.submit(self._process_organism, organism)
                for organism in self.organisms.values()
                ]
            for future in as_completed(futures):
                future.result()
        print("Finished building FASTA files!")