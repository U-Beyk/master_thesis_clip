'''
fasta_builder.py
================

Contains the class to build FASTA files.

author: U.B.
'''

from concurrent.futures import ProcessPoolExecutor, as_completed

from ..constants import ORGANISMS_TO_EXAMINE
from .rbp_site_generator import RbpSiteGenerator

class FastaBuilder:
    '''
    Class creating the FASTA files of specified organisms.

    Attributes
    ----------
    organisms: dict[str, str]
        Dictionary of all organisms to create FASTAs of.
    '''

    def __init__(self):
        '''
        Instantiates a FastaBuilder object.
        '''
        self.organisms: dict[str, str] = ORGANISMS_TO_EXAMINE

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
        Builds FASTA files of all organisms specified in the imported constant.
        The FASTA files are build separately each with its own worker.
        When finished a message is shown.
        '''
        with ProcessPoolExecutor(max_workers=7) as executor:
            futures = {
            executor.submit(self._process_organism, organism): name
            for name, organism in self.organisms.items()
        }
        for future in as_completed(futures):
            name = futures[future]
            try:
                future.result()
                print(f"Finished building FASTA for: {name}")
            except Exception as e:
                print(f"Error while processing {name}: {e}")
