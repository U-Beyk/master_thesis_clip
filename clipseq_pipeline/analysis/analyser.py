from typing import Callable

import pandas as pd

from .formatter import FormattedRbpDf
#
from .reporter import chi2_proteins_motifs
#

type DataFrameInitializer = Callable[[str, str], pd.DataFrame]
type DataFrameTransformer = Callable[[pd.DataFrame], pd.DataFrame]
type DataFrameFormatter = Callable[[pd.DataFrame], pd.DataFrame]
type DataFramePlotter= Callable[[pd.DataFrame, str], None]
type DataFrameReporter = Callable[[pd.DataFrame], None]

class AnalysisPipeline:
    def __init__(
        self,
        organism: str,
        algorithm: str,
        initializer: DataFrameInitializer,
        transformer: DataFrameTransformer,
        formatter: DataFrameFormatter,
        plotter: DataFramePlotter,
        #reporter: DataFrameReporter
    ):
        self.base_filepath = f"./data/plots/{organism}/{algorithm}"
        self.initializer = initializer
        self.transformer = transformer
        self.formatter = formatter
        self.plotter = plotter
        #self.reporter = reporter

    def _initialize_df(self, fasta_filepath: str, predictions_filepath: str) -> pd.DataFrame:
        """Create a DataFrame using the provided initializer."""
        return self.initializer(fasta_filepath, predictions_filepath)

    def _transform_df(self, df: pd.DataFrame) -> pd.DataFrame:
        """Apply the transformer function."""
        return self.transformer(df)
    
    def _format_dfs(self, dfs: pd.DataFrame) -> pd.DataFrame:
        return self.formatter(dfs)
    
    def _plot_dfs(self, dfs: pd.DataFrame) -> None:
        return self.plotter(dfs, self.base_filepath)
    
    #def _report_dfs(self, dfs: pd.DataFrame) -> None:
    #    return self.reporter(dfs)

    def run(self, fasta_filepath: str, predictions_filepath: str) -> None:
        """Full analysis pipeline."""
        df = self._initialize_df(fasta_filepath, predictions_filepath)
        df = self._transform_df(df)
        df = self._format_dfs(df)
        self._plot_dfs(df)
        #self._report_dfs(dfs)

        #Print for testing, remove later:
        #print(chi2_proteins_motifs(dfs))
        #