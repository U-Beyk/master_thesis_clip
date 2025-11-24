from typing import Callable

import pandas as pd

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
        reporter: DataFrameReporter
    ):
        self.organism = organism
        self.algorithm = algorithm
        self.initializer = initializer
        self.transformer = transformer
        self.formatter = formatter
        self.plotter = plotter
        self.reporter = reporter

    # TODO: Fix method.
    def _initialize_df(self, fasta_filepath: str, predictions_filepath: str) -> pd.DataFrame:
        fasta_original = f"./data/fasta_files/{self.organism}/{self.organism}_rbp_sites.fasta"
        predictions_original = f"./data/rna_predictions/{self.organism}/{self.organism}_prediction.csv"
        fasta_shuffled = f"./data/fasta_files/{self.organism}/{self.organism}_shuffled_rbp_sites.fasta"
        predictions_shuffled = f"./data/rna_predictions/{self.organism}/{self.organism}_shuffled_prediction.csv"
        return self.initializer(fasta_filepath, predictions_filepath)

    def _transform_df(self, df: pd.DataFrame) -> pd.DataFrame:
        return self.transformer(df)
    
    def _format_dfs(self, df: pd.DataFrame) -> pd.DataFrame:
        return self.formatter(df)
    
    def _plot_dfs(self, df: pd.DataFrame) -> None:
        base_filepath = f"./data/plots/{self.organism}/{self.algorithm}"
        return self.plotter(df, base_filepath)
    
    def _report_dfs(self, df: pd.DataFrame) -> None:
        base_filepath = f"./data/reports/{self.organism}/{self.algorithm}"
        return self.reporter(df, base_filepath)

    def run(self, fasta_filepath: str, predictions_filepath: str) -> None:
        df = self._initialize_df(fasta_filepath, predictions_filepath)
        df = self._transform_df(df)
        df = self._format_dfs(df)
        self._plot_dfs(df)
        self._report_dfs(df)