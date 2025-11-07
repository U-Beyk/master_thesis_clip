from .initializer import init_motifold_df, init_motices_df
from .transformer import trans_motifold_dfs, trans_motices_dfs
from .formatter import format_motifold_dfs, format_motices_dfs
from .plotter import plot_motifold_dfs

#
from.analyser import AnalysisPipeline

def run_analysis():
    analyser = AnalysisPipeline(
        "mus_musculus",
        "rnamotifold",
        init_motifold_df, 
        trans_motifold_dfs,
        format_motifold_dfs,
        plot_motifold_dfs
    )
    analyser.run(
        "./data/fasta_files/mus_musculus_rbp_sites.fasta",
        "./data/rna_predictions/mus_musculus_prediction.csv"
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