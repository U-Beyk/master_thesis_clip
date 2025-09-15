"""
data_mapper.py
~~~~~~~~~~~~~~

Contains all classes to map the CLIP data and sequences to specified transcripts.

Classes
-------
TranscriptProcessor
    Maps the CLIP data and sequence to their correpsonding transcript of a chromosome.
    Outputs the mapped data of a chromosome in a .json file.
PipelineRunner
    Pipeline for the mapping of the data for the specified organism.
SpeciesProcessor
    Maps all data for all chromosome of all specified organisms in the config.yaml file.

author: U.B.
"""

import yaml

from processing.rbp_site_generator import RbpSiteGenerator

class FastaBuilder:
    def __init__(self, config_file):
        self.organisms: dict[str, str] = self._extract_organisms(config_file)
    
    @staticmethod
    def _extract_organisms(config_file: str) -> dict[str, str]:
        with open(config_file, "r") as file:
            content = yaml.safe_load(file)
        return content["organisms"]
    
    def build_fasta_organisms(self) -> None:
        for _, organism in self.organisms.items():
            rbpsite_generator = RbpSiteGenerator(organism)
            rbpsite_generator.write_fasta()