from typing import Callable

import pandas as pd

type DataFrameInitializer = Callable[[str, str], pd.DataFrame]
type DataFrameTransformer = Callable[[pd.DataFrame], pd.DataFrame]
type DataFrameFormatter = Callable[[pd.DataFrame], pd.DataFrame]
type DataFramePlotter= Callable[[pd.DataFrame, str], None]
type DataFrameReporter = Callable[[pd.DataFrame, pd.DataFrame], None]
type ReportInspector = Callable[[str, str], None]

class AnalysisPipeline:
    def __init__(
        self,
        organism: str,
        algorithm: str,
        initializer: DataFrameInitializer,
        transformer: DataFrameTransformer,
        formatter: DataFrameFormatter,
        plotter: DataFramePlotter,
        reporter: DataFrameReporter,
        inspector: ReportInspector
    ):
        self.organism = organism
        self.algorithm = algorithm
        self.initializer = initializer
        self.transformer = transformer
        self.formatter = formatter
        self.plotter = plotter
        self.reporter = reporter
        self.inspector = inspector

    @property
    def fasta_original_path(self) -> str:
        return f"./data/fasta_files/{self.organism}/{self.organism}_rbp_sites.fasta"
    
    @property
    def predictions_original_path(self) -> str:
        return f"./data/rna_predictions/{self.organism}/{self.organism}_prediction.csv"
    
    @property
    def fasta_random_path(self) -> str:
        return f"./data/fasta_files/{self.organism}/{self.organism}_random_seq.fasta"
    
    @property
    def predictions_random_path(self) -> str:
        return f"./data/rna_predictions/{self.organism}/{self.organism}_random_preds.csv"
    
    @property
    def plot_directory_original(self) -> str:
        return f"./data/plots/{self.organism}/{self.algorithm}"
    
    @property
    def plot_directory_shuffled(self) -> str:
        return f"./data/plots/{self.organism}/nullmodel"
    
    @property
    def report_directory(self) -> str:
        return f"./data/reports/{self.organism}/{self.algorithm}"

    def _initialize_df(self, fasta_filepath: str, predictions_filepath: str) -> pd.DataFrame:
        return self.initializer(fasta_filepath, predictions_filepath)

    def _transform_df(self, df: pd.DataFrame) -> pd.DataFrame:
        return self.transformer(df)
    
    def _format_df(self, df: pd.DataFrame) -> pd.DataFrame:
        return self.formatter(df)
    
    def _plot_dfs(self, df: pd.DataFrame, directory_path: str) -> None:
        return self.plotter(df, directory_path)
    
    def _report_dfs(self, df: pd.DataFrame, df_null: pd.DataFrame) -> None:
        return self.reporter(df, df_null, self.report_directory)
    
    def _inspect_report(self) -> None:
        return self.inspector(self.organism, self.report_directory)

    def run(self) -> None:
        df = self._initialize_df(self.fasta_original_path, self.predictions_original_path)
        df_null = self._initialize_df(self.fasta_random_path, self.predictions_random_path)

        df = self._transform_df(df)
        df_null = self._transform_df(df_null)

        df = self._format_df(df)
        df_null = self._format_df(df_null)

        self._report_dfs(df, df_null)
        self._inspect_report()

        self._plot_dfs(df, self.plot_directory_original)
        self._plot_dfs(df_null, self.plot_directory_shuffled)