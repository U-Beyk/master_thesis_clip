'''
main.py
~~~~~~~

Main file which works as a pipeline for the processing, predicting and analysis of the CLIP data.

Functions
---------
main
    Main function running the processing, predicting and analysis of the data.
analyse_clip_data
    Analyses all CLIP data.
process_species
    Processes the data of all species by mappig the CLIP data and the sequences
    to their correspondinfg transcripts.
predict_species
    Predicts the structures of the sequences of the CLIP data.

author: U.B.
'''

from pre_processing_analysis.clip_analysis import ClipSpeciesAnalyzer
from processing.fasta_builder import FastaBuilder
from processing.rna_predicter import RNAPredicter
from prediction_analysis.species_analyzer import SpeciesAnalyzer

def main():
    '''
    Main function processing, predicting and analyzing the data.
    '''
    #process_species()
    #predict_species()
    analyze_predictions()

# TO-DO: Put these classes together into one???
def analyse_clip_data() -> None:
    '''
    Analyzes the CLIP data.
    '''
    analyzer = ClipSpeciesAnalyzer("config.yaml")
    analyzer.analyse_all_organisms()

def process_species() -> None:
    '''
    Processes the data of all species, specified in the config.yaml.
    '''
    fasta_builder = FastaBuilder("config.yaml")
    fasta_builder.build_fasta_organisms()

def predict_species() -> None:
    '''
    Predicts the structures and motifs of all organisms specified in the config.yaml.
    '''
    predicter = RNAPredicter("config.yaml")
    predicter.run_predictions()

def analyze_predictions() -> None:
    analyzer = SpeciesAnalyzer("config.yaml")
    analyzer.analyze_all_organisms()

if __name__ == "__main__":
    main()