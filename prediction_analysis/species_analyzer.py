import yaml
from concurrent.futures import ProcessPoolExecutor, as_completed

from prediction_analysis.rna_predictions_visualizer import CLIPPredictionVisualizer, AnalyzerConfig

class SpeciesAnalyzer:

    def __init__(self, config_file: str):
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
    def _analyze_organism(organism: str) -> None:
        config = AnalyzerConfig(
        f"./data/fasta_files/{organism}_rbp_sites.fasta",
        f"./data/rna_predictions/{organism}_prediction.csv"
        )
        visualizer = CLIPPredictionVisualizer(
            config,
            organism
        )
        visualizer.visualize_clip_data()

    def analyze_all_organisms(self) -> None:
        print("Analyzing predictions!")
        with ProcessPoolExecutor(max_workers=7) as executor:
            futures = [
                executor.submit(self._analyze_organism, organism)
                for organism in self.organisms.values()
                ]
            for future in as_completed(futures):
                future.result()
        print("Finished the analysis!")