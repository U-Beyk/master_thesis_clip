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

from processing.fasta_builder import FastaBuilder
from pre_processing_analysis.clip_analysis import ClipSpeciesAnalyzer

def main():
    """
    Main function processing, predicting and analyzing the data. 
    """
    process_species()
    #analyse_clip_data()

def analyse_clip_data() -> None:
    analyzer = ClipSpeciesAnalyzer("config.yaml")
    analyzer.analyse_all_organisms()

def process_species() -> None:
    """
    Processes the data of all species, specified in the config.yaml.
    """
    fasta_builder = FastaBuilder("config.yaml")
    fasta_builder.build_fasta_organisms()

if __name__ == "__main__":
    main()