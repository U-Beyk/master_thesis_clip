"""
main.py

Main file which works as a pipeline for the processing , predicting and analysis of the CLIP data.

Functions
---------
main
    Main function running the processing, predicting and analysis of the data.
process_species
    Processes the data of all species by mappig the CLIP data and the sequences
    to their correspondinfg transcripts.

author: U.B.
"""

from processing.data_mapper import SpeciesProcessor

def main():
    """
    Main function processing, predicting and analyzing the data. 
    """
    process_species()

def process_species() -> None:
    """
    Processes the data of all species, specified in the config.yaml.
    """
    processor = SpeciesProcessor("config.yaml")
    processor.process_species()

if __name__ == "__main__":
    main()