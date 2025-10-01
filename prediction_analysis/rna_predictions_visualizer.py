'''
Classes using the data to visualize the results of the analysis.

author: U.B.
'''

import pandas as pd
from dataclasses import dataclass, field

from prediction_analysis.analyzer import AnalyzerConfig, Analyzer
from prediction_analysis.records_assembler import CLIPmotiFoldRecordsAssembler, CLIPmotiCesRecordsAssembler
from prediction_analysis.visualization import HeatMap

#TO-DO: Add RNAmotiCes prediction visualization and heatmap
@dataclass
class PlotConfig:
    '''Class for the settings of the plots.'''
    organism: str
    folder: str
    extra_info: str
    base_plot_path: str = field(init=False, default="./plots")

    def shape_plot_path(self, plot_type: str) -> str:
        '''Creates the path for the shape abstraction files.'''
        return f"{self.base_plot_path}/shape_abstraction/{self.organism}/{self.folder}/{self.organism}_{plot_type}.jpg"
    
class MotifoldPlotter:
    '''Class plotting the RNAmotiFOld results.'''
    def __init__(self, assembly: CLIPmotiFoldRecordsAssembler, config: PlotConfig) -> None:
        '''Initializes a plotter object for the RNAmotiFold assemblies.'''
        self.assembly = assembly
        self.config = config

    def barchart(self) -> None:
        '''Creates and saves a barchart of the RNAmotiFold assembly.'''
        self.assembly.visualize_as_barchart(
            f"{self.config.organism}: Motif occurrences per sequence\n{self.config.extra_info}",
            self.config.shape_plot_path("barchart")
        )

    def violinplot(self) -> None:
        '''Creates and saves a violinplot of the RNAmotiFold assembly.'''
        self.assembly.visualize_as_violinplot(
            f"{self.config.organism}: Distance of mfe to lowest mfe\n{self.config.extra_info}",
            self.config.shape_plot_path("violinplot")
        )

class CLIPPredictionVisualizer:
    '''Class visualizing the CRISPR RNA predictions.'''
    def __init__(self, config: AnalyzerConfig, organism: str) -> None:
        '''Initializes a CRISPRRNAPredictionVisualizer object.'''
        self.analyzer = Analyzer(config)
        self.organism = organism

    def _visualize_all_data_mfe_below_zero(self) -> None:
        '''Visualize the data for all RNAs for mfe<0'''
        config = PlotConfig(self.organism, "all_data", "all CLIP data, mfe < 0 ")
        MotifoldPlotter(self.analyzer.rna_motifold, config).barchart()
        MotifoldPlotter(self.analyzer.rna_motifold, config).violinplot()

    # TO-DO: Check if filtering runs fine.
    def _visualize_unique_proteins(self) -> None:
        for rbp in self.analyzer.rna_motifold.unique_rbps:
            motifold = self.analyzer.filter_motifold_by_rbp(rbp)
            config = PlotConfig(self.organism, f"{rbp.lower()}", f"Clip data of {rbp}, mfe < 0")
            MotifoldPlotter(motifold, config).barchart()
            MotifoldPlotter(motifold, config).violinplot()

    def visualize_clip_data(self) -> None:
        '''Visualize the data for all RNAs for mfe<0 and within the mfe range.'''
        self._visualize_all_data_mfe_below_zero()
        self._visualize_unique_proteins()