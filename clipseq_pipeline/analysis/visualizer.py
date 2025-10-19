from dataclasses import dataclass, field

from .statistics import CLIPmotiFoldStatistics, CLIPmotiCesStatistics
from .visualization import BarChart, ViolinPlot, Histogram, HeatMap

@dataclass
class PlotConfig:
    '''Class for the settings of the plots.'''
    organism: str
    folder: str
    extra_info: str
    base_plot_path: str = field(init=False, default="./data/plots")

    def shape_plot_path(self, plot_type: str) -> str:
        '''Creates the path for the shape abstraction files.'''
        return f"{self.base_plot_path}/{self.organism}/shape_abstraction/{self.folder}/{self.organism}_{plot_type}.jpg"

class CLIPmotiFoldVisualizer:
    def __init__(self, stats: CLIPmotiFoldStatistics, config:PlotConfig):
        self.stats = stats
        self.config = config

    def visualize_as_barchart(self) -> None:
        barchart = BarChart(
            f"{self.config.organism}: Motif occurrences per sequence\n{self.config.extra_info}",
            self.stats.clip_overlap_motifs_count,
            self.stats.sequence_number
        )
        barchart.save_plot(self.config.shape_plot_path("barchart"))

    def visualize_as_violinplot(self) -> None:
        violinplot = ViolinPlot(
            f"{self.config.organism}: Distance of mfe to lowest mfe\n{self.config.extra_info}",
            self.stats.mfe_clip_motif_distances,
            self.stats.prediction_number
        )
        violinplot.save_plot(self.config.shape_plot_path("violinplot"))

# TODO: Adjust class, especially histogram method
class CLIPmotiCesVisualizer(CLIPmotiFoldVisualizer):
    def __init__(self, stats: CLIPmotiCesStatistics, config: PlotConfig):
        self.stats = stats
        self.config

    def visualize_as_histogram(self, description: str, filepath: str) -> None:
        histogram = Histogram(
            description,
            self.stats.motifs_count,
            len(self.stats.rna_sequences)
        )
        histogram.save_plot(filepath)