from typing import Callable

import pandas as pd

from .transformer import FilteredRbpDf
from .formatter import FormattedRbpDf
#
from .reporter import chi2_proteins_motifs
#

type DataFrameInitializer = Callable[[str, str], pd.DataFrame]
type DataFrameTransformer = Callable[[pd.DataFrame], list[FilteredRbpDf]]
type DataFrameFormatter = Callable[[list[FilteredRbpDf]], list[FormattedRbpDf]]
type DataFramePlotter= Callable[[list[FormattedRbpDf], str]]

class AnalysisPipeline:
    def __init__(
        self,
        organism: str,
        algorithm: str,
        initializer: DataFrameInitializer,
        transformer: DataFrameTransformer,
        formatter: DataFrameFormatter,
        plotter: DataFramePlotter
    ):
        self.base_filepath = f"./data/plots/{organism}/{algorithm}"
        self.initializer = initializer
        self.transformer = transformer
        self.formatter = formatter
        self.plotter = plotter

    def _initialize_df(self, fasta_filepath: str, predictions_filepath: str) -> pd.DataFrame:
        """Create a DataFrame using the provided initializer."""
        return self.initializer(fasta_filepath, predictions_filepath)

    def _transform_df(self, df: pd.DataFrame) -> list[FilteredRbpDf]:
        """Apply the transformer function."""
        return self.transformer(df)
    
    def _format_dfs(self, dfs: list[FilteredRbpDf]) -> list[FormattedRbpDf]:
        return self.formatter(dfs)
    
    def _plot_dfs(self, dfs: list[FormattedRbpDf]) -> None:
        return self.plotter(dfs, self.base_filepath)

    def run(self, fasta_filepath: str, predictions_filepath: str) -> None:
        """Full analysis pipeline."""
        df = self._initialize_df(fasta_filepath, predictions_filepath)
        dfs = self._transform_df(df)
        dfs = self._format_dfs(dfs)
        #self._plot_dfs(dfs)
        #
        print(chi2_proteins_motifs(dfs))
        #