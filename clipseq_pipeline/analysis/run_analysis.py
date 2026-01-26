from .species_analyser import SpeciesAnalyser


def run_analysis():
    analyser = SpeciesAnalyser()
    analyser.analyse_organisms()