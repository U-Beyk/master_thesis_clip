'''
Contains all classes for the visualization of the RNA prediction data.

author: U.B.
'''

from abc import ABC, abstractmethod
from dataclasses import dataclass
import os
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from matplotlib.axes import Axes
from matplotlib.ticker import FixedLocator

# TODO: Change parent class to Plot and filename to plot.py.
# TODO: Refactor classes.
@dataclass
class Visualization(ABC):
    '''Abstract base class for motif analysis visualizations.'''
    title: str

    @abstractmethod
    def _draw(self, ax: Axes) -> None:
        '''Abstract method for drawing the visualization.'''
        pass

    def save_plot(self, filepath: str) -> None:
        '''Saves the generated plot.'''
        os.makedirs(os.path.dirname(filepath), exist_ok=True)
        fig, ax = plt.subplots(figsize=(12, 10))
        self._draw(ax)
        fig.tight_layout()
        fig.savefig(f"{filepath}", dpi=300)
        plt.close(fig)

@dataclass
class BarChart(Visualization):
    '''Bar chart representation of motif analysis results.'''
    motif_numbers: pd.Series
    max_value: int

    def _draw(self, ax: Axes) -> None:
        '''Draws bar chart visualization.'''
        sns.barplot(
            x=self.motif_numbers.index,
            y=self.motif_numbers.values,
            ax=ax,
            color="blue"
        )
        for i, value in enumerate(self.motif_numbers.values):
            ax.text(
                i, value, str(int(value)),
                ha="center", va="bottom", fontsize=16
            )
        self._style_axes(ax)

    def _add_entry_count(self, ax: Axes) -> None:
        '''Adds the total number of entries to the plot at the top right corner.'''
        ax.text(
            0.95, 0.95, f"Number of sequences={self.max_value}", 
            transform=ax.transAxes, 
            fontsize=20,
            verticalalignment='top', 
            horizontalalignment='right',
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.7)
        )

    def _style_axes(self, ax: Axes) -> None:
        '''Styles bar chart axes.'''
        max_height = self.motif_numbers.max()
        ax.set_ylim(0, max_height * 1.15)
        ax.set_title(self.title, pad=20)
        ax.set_xlabel("Motifs", fontsize=20)
        ax.set_ylabel("Count", fontsize=18)
        ax.tick_params(axis="x", labelrotation=45, labelsize=22)
        ax.tick_params(axis="y", labelsize=20)
        plt.setp(ax.get_xticklabels(), ha="right")
        self._add_entry_count(ax)

# TODO: Refactor Histogram, so that it takes a pd.Series
@dataclass
class Histogram(BarChart):
    '''Histogram representation, inherits functionality from BarChart.'''
    def _draw(self, ax: Axes) -> None:
        '''Draws histogram visualization.'''
        ax.bar(self.motifs, self.numbers, color="blue", width=0.8)
        # Add labels
        for bar in ax.patches:
            yval = bar.get_height()
            ax.text(bar.get_x() + bar.get_width() / 2, yval, str(int(yval)),
                    ha="center", va="bottom", fontsize=16)

        self._style_axes(ax)
    
    def _style_axes(self, ax: Axes) -> None:
        '''Styles histogram axes. Scaling of the bars is adjusted to the highest bar.'''
        max_height = max(self.numbers)
        ax.set_ylim(0, max_height * 1.15)
        ax.set_title(self.title, pad=20)
        ax.set_xlabel("Sequences", fontsize=20)
        ax.tick_params(axis="x", labelrotation=45)
        ax.tick_params(axis="y", labelsize=20)
        ax.set_xticks(range(len(self.motifs)))
        ax.set_xticklabels(self.motifs, ha="right", fontsize=22)
        self._add_entry_count(ax)

@dataclass
class ViolinPlot(Visualization):
    '''Base class for violin and box plots.'''
    energies_motifs: pd.DataFrame
    max_value: int

    def _style_axes(self, ax: Axes) -> None:
        '''Configures plot styling.'''
        ax.set_title(self.title)
        ax.set_xlabel("Motifs", fontsize=20)
        ax.set_ylabel("Distance to mfe [kcal/mol]", fontsize=18)
        ax.tick_params(axis="x", labelrotation=60)
        ax.tick_params(axis="y", labelsize=20)
        ax.grid(True, axis="y", linestyle="--", linewidth=0.5, alpha=0.7)
        self._add_counts_to_labels(ax)
        self._add_entry_count(ax)

    def _add_counts_to_labels(self, ax: Axes) -> None:
        '''Adds the count of each motif to the x-axis labels.'''
        ax.figure.canvas.draw()
        counts = self.energies_motifs['Motifs'].value_counts()
        xticks = ax.get_xticks()
        xticklabels = [f"{label.get_text()}\n(n={counts.get(label.get_text(), 0)})" for label in ax.get_xticklabels()]
        ax.xaxis.set_major_locator(FixedLocator(xticks))
        ax.set_xticklabels(xticklabels, fontsize=22)

    def _add_entry_count(self, ax: Axes) -> None:
        '''Adds the total number of entries to the plot.'''
        ax.text(
            0.95, 0.95, f"Number of predictions={self.max_value}", 
            transform=ax.transAxes, 
            fontsize=20,
            verticalalignment='top', 
            horizontalalignment='right',
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.7)
        )

    def _draw(self, ax: Axes) -> None:
        '''Base draw method.'''
        sns.violinplot(
            x="Motifs", y="Distance to mfe", hue="Motifs", data=self.energies_motifs,
            inner="box", density_norm="width", palette="Set2", legend=False, ax=ax, cut=0)
        self._style_axes(ax)

# TODO: Refactor Heatmap class.
@dataclass
class HeatMap(Visualization):
    '''Heatmap representation of motif median values.'''
    type_motif_medians: pd.DataFrame

    def _reorder_motifs(self) -> None:
        '''Reorders motifs based on non-null counts.'''
        motif_order = self.type_motif_medians.notna().sum(axis=1).sort_values(ascending=False).index
        motif_order = motif_order.drop("No motif")
        motif_order = ["No motif"] + list(motif_order)
        self.type_motif_medians = self.type_motif_medians.loc[motif_order]

    def _draw(self, ax: Axes) -> None:
        '''Creates heatmap visualization.'''
        self._reorder_motifs()
        mask = self.type_motif_medians.isna()
        heatmap = sns.heatmap(self.type_motif_medians, annot=True, cbar=True, fmt=".2f", linewidths=0.5,
                    linecolor="gray", mask=mask, ax=ax, annot_kws={"fontsize": 13})
        self._style_axes(ax, heatmap)

    def _style_axes(self, ax: Axes, heatmap: Axes) -> None:
        '''Styles heatmap axes.'''
        ax.set_title(self.title)
        ax.tick_params(axis='x', labelsize=18)
        ax.tick_params(axis='y', labelsize=20, labelrotation=0)
        ax.set_ylabel("")
        cbar = heatmap.collections[0].colorbar
        cbar.ax.tick_params(labelsize=14)