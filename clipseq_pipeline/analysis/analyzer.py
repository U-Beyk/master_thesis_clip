'''
Contains all classes with the dataframes to be analyzed and its methods to filter necessary information.

author: U.B.
'''

import copy
from dataclasses import dataclass

from .assembler import CLIPmotiFoldAssembler
from .statistics import CLIPmotiFoldStatistics
from .visualizer import PlotConfig, CLIPmotiFoldVisualizer, CLIPmotiCesVisualizer

# TODO: Classes have to be refactored.
@dataclass
class AnalyzerConfig:
    '''Dataclass containing the information of the fasta file and the RNAmotiFold and RNAHeliCes/RNAmotiCes prediction files.'''
    fasta_path: str
    motifold_csv_path: str

class Analyzer:
    '''Generic class for the different RNA assemblies.'''

    def __init__(self, config: AnalyzerConfig, organism: str) -> None:
        """Initializes an Analyzer object with an RNAmotiFold assembly."""
        self.organism = organism
        self.rna_motifold = CLIPmotiFoldAssembler(config.fasta_path, config.motifold_csv_path)
        self.motifold_stats = CLIPmotiFoldStatistics(self.rna_motifold.rna_sequences)

    def _filter_motifold_mfe_range(self) -> CLIPmotiFoldAssembler:
        """Filters RNAmotiFold predictions within an mfe range (10% above lowest mfe)."""
        threshold = self.motifold_stats.lowest_mfe_value * 0.1
        return self.rna_motifold.filter_records(mfe_threshold=threshold)
    
    def _filter_motifold_lowest_mfe(self) -> CLIPmotiFoldAssembler:
        return self.rna_motifold.filter_lowest_mfe_predictions()

    def _filter_motifold_by_rbp(self, rbp: str) -> CLIPmotiFoldAssembler:
        """Filters RNAmotiFold records by a specific RBP, 
        optionally restricting to an mfe range."""
        return self.rna_motifold.filter_records(lambda r: r.rbp_name == rbp)
    
    def _filter_motifold_by_rbp_lowest_mfe(self, rbp: str) -> CLIPmotiFoldAssembler:
        filtered = self.rna_motifold.filter_lowest_mfe_predictions()
        return filtered.filter_records(lambda r: r.rbp_name == rbp)
    
    def _filter_motifold_by_feature_mrna(self):
        return self.rna_motifold.filter_records(lambda r: "mRNA" in r.feature_types)
    
    def _filter_motifold_by_feature_ncrna(self):
        return self.rna_motifold.filter_records(lambda r: "ncRNA_gene" in r.feature_types)
    
    def _visualize_all_data_mfe_below_zero(self) -> None:
        """Visualize all RNAs for mfe < 0 using barchart and violinplot."""
        config = PlotConfig(
            self.organism,
            "all_data",
            "All CLIP data, mfe < 0"
        )
        plotter = CLIPmotiFoldVisualizer(self.motifold_stats, config)
        plotter.visualize_as_barchart()
        plotter.visualize_as_violinplot()

    def _visualize_all_data_mfe_range(self) -> None:
        """Visualize all RNAs for mfe < 0 using barchart and violinplot."""
        filtered_assembler = self._filter_motifold_mfe_range()
        filtered_stats = CLIPmotiFoldStatistics(filtered_assembler.rna_sequences)
        config = PlotConfig(
            self.organism,
            "all_data",
            "All CLIP data, mfe range"
        )
        plotter = CLIPmotiFoldVisualizer(filtered_stats, config)
        plotter.visualize_as_barchart()
        plotter.visualize_as_violinplot()

    def _visualize_all_data_lowest_mfe(self) -> None:
        """Visualize all RNAs for mfe < 0 using barchart and violinplot."""
        filtered_assembler = self._filter_motifold_lowest_mfe()
        filtered_stats = CLIPmotiFoldStatistics(filtered_assembler.rna_sequences)
        config = PlotConfig(
            self.organism,
            "all_data",
            "All CLIP data, lowest mfe"
        )
        plotter = CLIPmotiFoldVisualizer(filtered_stats, config)
        plotter.visualize_as_barchart()
        plotter.visualize_as_violinplot()

    def _visualize_all_by_mrna(self) -> None:
        """Visualize all RNAs for mfe < 0 using barchart and violinplot."""
        filtered_assembler = self._filter_motifold_by_feature_mrna()
        filtered_stats = CLIPmotiFoldStatistics(filtered_assembler.rna_sequences)
        config = PlotConfig(
            self.organism,
            "all_data",
            "CLIP feature mRNA, mfe < 0"
        )
        plotter = CLIPmotiFoldVisualizer(filtered_stats, config)
        plotter.visualize_as_barchart()
        plotter.visualize_as_violinplot()

    def _visualize_all_by_ncrna(self) -> None:
        """Visualize all RNAs for mfe < 0 using barchart and violinplot."""
        filtered_assembler = self._filter_motifold_by_feature_ncrna()
        filtered_stats = CLIPmotiFoldStatistics(filtered_assembler.rna_sequences)
        config = PlotConfig(
            self.organism,
            "all_data",
            "CLIP feature ncRNA, mfe < 0"
        )
        plotter = CLIPmotiFoldVisualizer(filtered_stats, config)
        plotter.visualize_as_barchart()
        plotter.visualize_as_violinplot()

    def _visualize_unique_proteins(self) -> None:
        """Visualize each unique RBP individually."""
        for rbp in self.motifold_stats.unique_rbps:
            filtered_assembler = self._filter_motifold_by_rbp(rbp)
            filtered_stats = CLIPmotiFoldStatistics(filtered_assembler.rna_sequences)
            config = PlotConfig(self.organism, f"{rbp.lower()}", f"CLIP data of {rbp}, mfe < 0")
            plotter = CLIPmotiFoldVisualizer(filtered_stats, config)
            plotter.visualize_as_barchart()
            plotter.visualize_as_violinplot()

    def _visualize_unique_proteins_mfe_range(self) -> None:
        """Visualize each unique RBP individually."""
        for rbp in self.motifold_stats.unique_rbps:
            filtered_assembler = self._filter_motifold_by_rbp(rbp, mfe_range=True)
            filtered_stats = CLIPmotiFoldStatistics(filtered_assembler.rna_sequences)
            config = PlotConfig(self.organism, f"{rbp.lower()}", f"CLIP data of {rbp}, mfe range")
            plotter = CLIPmotiFoldVisualizer(filtered_stats, config)
            plotter.visualize_as_barchart()
            plotter.visualize_as_violinplot()

    def _visualize_unique_proteins_lowest_mfe(self) -> None:
        """Visualize each unique RBP individually."""
        for rbp in self.motifold_stats.unique_rbps:
            filtered_assembler = self._filter_motifold_by_rbp_lowest_mfe(rbp)
            filtered_stats = CLIPmotiFoldStatistics(filtered_assembler.rna_sequences)
            config = PlotConfig(self.organism, f"{rbp.lower()}", f"CLIP data of {rbp}, lowest mfe")
            plotter = CLIPmotiFoldVisualizer(filtered_stats, config)
            plotter.visualize_as_barchart()
            plotter.visualize_as_violinplot()

    def _visualize_unique_proteins_by_mrna(self) -> None:
        """Visualize each unique RBP individually."""
        for rbp in self.motifold_stats.unique_rbps:
            filtered_assembler = self._filter_motifold_by_feature_mrna().filter_records(lambda r: r.rbp_name == rbp)
            filtered_stats = CLIPmotiFoldStatistics(filtered_assembler.rna_sequences)
            config = PlotConfig(self.organism, f"{rbp.lower()}", f"CLIP data of {rbp}, mfe < 0")
            plotter = CLIPmotiFoldVisualizer(filtered_stats, config)
            plotter.visualize_as_barchart()
            plotter.visualize_as_violinplot()

    def _visualize_unique_proteins_by_ncrna(self) -> None:
        """Visualize each unique RBP individually."""
        for rbp in self.motifold_stats.unique_rbps:
            filtered_assembler = self._filter_motifold_by_feature_ncrna().filter_records(lambda r: r.rbp_name == rbp)
            filtered_stats = CLIPmotiFoldStatistics(filtered_assembler.rna_sequences)
            config = PlotConfig(self.organism, f"{rbp.lower()}", f"CLIP data of {rbp}, mfe < 0")
            plotter = CLIPmotiFoldVisualizer(filtered_stats, config)
            plotter.visualize_as_barchart()
            plotter.visualize_as_violinplot()

    def visualize_clip_data(self) -> None:
        """Visualize all data and individual RBPs."""
        self._visualize_all_by_ncrna()
        self._visualize_unique_proteins_by_ncrna()

    def __repr__(self) -> str:
        return (
            f"{self.__class__.__name__}("
            f"rna_motifold_assembly={self.rna_motifold})"
        )