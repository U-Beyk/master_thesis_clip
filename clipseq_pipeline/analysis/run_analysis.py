from .species_analyzer import SpeciesAnalyzer

def run_analysis():
    analyser = SpeciesAnalyzer()
    analyser.analyze_all_organisms()