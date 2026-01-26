from.analyser import AnalysisPipeline
from ..constants import ORGANISMS_TO_EXAMINE
from .formatter import format_motifold_df, format_motices_df
from .initializer import init_motifold_df, init_motices_df
from .inspector import inspect_motifold_df
from .plotter import plot_motifold_df
from .reporter import report_motifold_df
from .transformer import trans_motifold_df, trans_motices_df

# TODO: Add multiprocessing for each organism!
class SpeciesAnalyser:

    def __init__(self):
        self.organisms: dict[str, str] = ORGANISMS_TO_EXAMINE

    def analyse_organisms(self) -> None:
        for name, organism  in self.organisms.items():
            self._analyse_rnamotifold_organism(name, organism)
    
    def _analyse_rnamotifold_organism(self, name: str, organism: str) -> None:
        analyser = AnalysisPipeline(
            organism,
            "rnamotifold",
            init_motifold_df, 
            trans_motifold_df,
            format_motifold_df,
            plot_motifold_df,
            report_motifold_df,
            inspect_motifold_df
        )
        analyser.run()
        print(f"Analysed RNAmotiFold results of {name}!")