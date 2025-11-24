from .initializer import init_motifold_df, init_motices_df
from .transformer import trans_motifold_df, trans_motices_df
from .formatter import format_motifold_df, format_motices_df
from .plotter import plot_motifold_df
from .reporter import report_motifold_df

#
from.analyser import AnalysisPipeline

def run_analysis():
    analyser = AnalysisPipeline(
        "drosophila_melanogaster",
        "rnamotifold",
        init_motifold_df, 
        trans_motifold_df,
        format_motifold_df,
        plot_motifold_df,
        report_motifold_df
    )
    analyser.run(
        "./data/fasta_files/drosophila_melanogaster_mono_nt_shuffled_rbp_sites.fasta",
        "./data/rna_predictions/drosophila_melanogaster_mono_nt_shuffled_prediction.csv"
    )
    #df = init_motifold_df(
    #    "./data/fasta_files/danio_rerio_rbp_sites.fasta",
    #    "./data/rna_predictions/danio_rerio_prediction.csv"
    #)
    #formatter = BarplotFormatter()
    #print(formatter.format(df))
    #formatted_df = formatter.format(df)
    #plotter = Barplot()
    #plotter.plot(formatted_df, "./data/plots.png")
    #dfs = transform_dfs(df)
    #print(dfs)