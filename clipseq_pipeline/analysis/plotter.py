from abc import ABC, abstractmethod
from dataclasses import dataclass
import os
from typing import Callable

import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator
import pandas as pd
import seaborn as sns

type PlotFn = Callable[[pd.DataFrame], None]

@dataclass
class RnaPlotter:
    plot_name: str
    plot_fn: PlotFn

class PlotStrategy(ABC):
    """Abstract base class for plot strategies."""

    @abstractmethod
    def plot(self, df: pd.DataFrame, filepath: str) -> None:
        """Generate a plot based on the given DataFrame and save it to a file."""
        pass

    def _save_plot(self, filepath: str) -> None:
        """Ensure directory exists and save the current matplotlib figure."""
        os.makedirs(os.path.dirname(filepath), exist_ok=True)
        plt.savefig(filepath, dpi=300, bbox_inches="tight")
        plt.close()


class Barplot(PlotStrategy):
    """Concrete strategy for generating bar plots."""

    def plot(self, df: pd.DataFrame, filepath: str) -> None:
        fig, ax = plt.subplots(figsize=(8, 5))
        sns.barplot(data=df, x="Motifs", y="Count")
        self._add_category_numbers(df, ax)
        self._save_plot(filepath)

    def _add_category_numbers(self, df: pd.DataFrame, ax: plt.Axes) -> None:
        for i, value in enumerate(df["Count"]):
            ax.text(
                i, value, str(int(value)),
                ha="center", va="bottom"
            )


class Violinplot(PlotStrategy):
    """Concrete strategy for generating violin plots."""

    def plot(self, df: pd.DataFrame, filepath: str) -> None:
        fig, ax = plt.subplots(figsize=(8, 5))
        self._draw_violinplot(df, ax)
        self._format_axes(df, ax)
        self._save_plot(filepath)

    def _draw_violinplot(self, df: pd.DataFrame, ax: plt.Axes) -> None:
        """Draw the main violin plot."""
        sns.violinplot(
            data=df,
            x="Motifs",
            y="Mfe distances [kcal/mol]",
            hue="Motifs",
            cut=0,
            ax=ax
        )

    def _format_axes(self,df: pd.DataFrame, ax: plt.Axes) -> None:
        """Format grid, labels, and ticks for clarity."""
        ax.grid(True, axis="y", linestyle="--", linewidth=0.5, alpha=0.5)
        ax.tick_params(axis="x", labelrotation=60)
        ax.figure.canvas.draw()

        counts = df["Motifs"].value_counts()
        categories = sorted(df["Motifs"].unique())
        xticklabels = [
            f"{category}\n(n={counts.get(category, 0)})"
            for category in categories
        ]
        xticks = range(len(categories))

        ax.xaxis.set_major_locator(FixedLocator(xticks))
        ax.set_xticklabels(xticklabels, fontsize=12)

# TODO: Implement class.
class Histogram(PlotStrategy):
    def plot(self, df: pd.DataFrame, filepath: str):
        raise NotImplementedError("Histogram Class and plot method are not implemented yet.")

# TODO: Refactor code.
def apply_plotters(
    formatted_df: pd.DataFrame,
    plotters: dict[str, RnaPlotter],
    base_filepath: str
) -> None:
    for (format_name, filter_name, rbp_name), row in formatted_df.iterrows():
        if format_name not in plotters:
            continue

        plotter = plotters[format_name]
        plot_fn = plotter.plot_fn

        df: pd.DataFrame = row["data"]
        if df.empty:
            continue

        output_dir = os.path.join(base_filepath, filter_name, rbp_name)
        os.makedirs(output_dir, exist_ok=True)
        output_path = os.path.join(output_dir, f"{format_name}_{plotter.plot_name}.png")

        plot_fn(df, output_path)

RNAMOTIFOLD_PLOTS = {
    "motif_frequency": RnaPlotter("barplot", Barplot().plot),
    "motif_distances": RnaPlotter("violinplot", Violinplot().plot)
}

def plot_motifold_df(formatted_df: pd.DataFrame, base_filepath: str) -> None:
    return apply_plotters(formatted_df, RNAMOTIFOLD_PLOTS, base_filepath)

RNAMOTICES_PLOTS = {
    "potential_motifs": RnaPlotter("hisotogram", Histogram().plot)
}

def plot_motices_df(formatted_df: pd.DataFrame, base_filepath: str) -> None:
    return apply_plotters(formatted_df, RNAMOTICES_PLOTS, base_filepath)