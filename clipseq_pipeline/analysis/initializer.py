import re

from Bio import SeqIO
import pandas as pd

from ..constants import CLIP_NT_WINDOW, AMBIGUOUS_MOTIFS

class DataLoader:
    def __init__(self, fasta_filepath: str, predictions_filepath: str):
        self.fasta_filepath = fasta_filepath
        self.predictions_filepath = predictions_filepath

    def build_dataframe(self) -> pd.DataFrame:
        sequences_df = self._generate_fasta_df()
        predictions_df = self._generate_prediction_df()
        return pd.merge(predictions_df, sequences_df, on="ID", how="inner")

    def _generate_fasta_df(self) -> pd.DataFrame:
        records = []
        for record in SeqIO.parse(self.fasta_filepath, "fasta"):
            records.append({"ID": record.id, "sequence": str(record.seq)})
        return pd.DataFrame(records)

    def _generate_prediction_df(self) -> pd.DataFrame:
        return pd.read_csv(self.predictions_filepath, sep="\t")

def parse_id_column(df: pd.DataFrame) -> pd.DataFrame:
    # Splits the ID column by "|" first
    split_id_df = df["ID"].str.split("|", expand=True)
    # Defines column names and their positions in the split
    columns_map = {
        "sequence_id": 0,
        "clip_range": 1,
        "rbp_name": 2,
        "seq_range": 3,
        "feature_types": 4,
    }
    # Extracts each column
    for col_name, index in columns_map.items():
        df[col_name] = split_id_df[index].str.split(":", expand=True).iloc[:, -1]
    # Splits ranges into start and end
    df[["clip_start", "clip_end"]] = df["clip_range"].str.split("-", expand=True).astype(int)
    df[["sequence_start", "sequence_end"]] = df["seq_range"].str.split("-", expand=True).astype(int)
    # Drops the now redundant columns
    df = df.drop(columns=["ID", "clip_range", "seq_range"])
    return df

def reorder_columns(df: pd.DataFrame) -> pd.DataFrame:
    new_cols = [
        "sequence_id", "clip_start", "clip_end","rbp_name", 
        "sequence_start", "sequence_end", "feature_types"
    ]
    df = df[new_cols + [col for col in df.columns if col not in new_cols]]
    return df

def adjust_mfe_to_kcal(df: pd.DataFrame) -> pd.DataFrame:
    df["mfe"] = df["mfe"] / 100
    return df

def add_relative_start_end(df: pd.DataFrame) -> pd.DataFrame:
    df["relative_start"] = df["clip_start"] - df["sequence_start"]
    df["relative_end"] = df["clip_end"] - df["sequence_start"]
    return df

def add_distance_to_mfe(df: pd.DataFrame) -> pd.DataFrame:
    df["mfe_distance"] = df.groupby("sequence_id")["mfe"].transform(lambda x: x - x.min())
    return df

def rename_columns(df: pd.DataFrame) -> pd.DataFrame:
    return df.rename(columns={
        "Class": "prediction_class",
        "motBracket": "mot_bracket"
    })

def add_clip_window(df: pd.DataFrame, target_length: int) -> pd.DataFrame:

    def extract_window(row):
        start, end = row["relative_start"], row["relative_end"]
        clip_len = end - start
        extend = max(0, target_length - clip_len)
        extend_left = extend // 2
        extend_right = extend - extend_left
        # Compute new start/end with boundaries
        seq_start = max(0, start - extend_left)
        seq_end = min(len(row["mot_bracket"]), end + extend_right)
        # If we couldn't extend fully on one side, push remaining to the other
        missing = target_length - (seq_end - seq_start)
        if missing > 0:
            seq_start = max(0, seq_start - missing)
            seq_end = min(len(row["mot_bracket"]), seq_end + missing)

        return row["mot_bracket"][seq_start:seq_end]

    df["mot_bracket_window"] = df.apply(extract_window, axis=1)
    return df

def replace_ambiguous_motifs(df: pd.DataFrame, ambiguous_motifs: dict[str, str]) -> pd.DataFrame:
    mask = df["prediction_class"].notna()
    df.loc[mask, "prediction_class"] = df.loc[mask, "prediction_class"].replace(ambiguous_motifs, regex=True)
    df["motifs_clip_window"] = df["motifs_clip_window"].replace(ambiguous_motifs, regex=True)
    return df

def fill_nan_motifs(df: pd.DataFrame) -> pd.DataFrame:
    df["prediction_class"] = df["prediction_class"].fillna("")
    return df

def add_motifs_clip_window(df: pd.DataFrame) -> pd.DataFrame:
    df["motifs_clip_window"] = df["mot_bracket_window"].apply(lambda x: "".join(set(re.findall(r"[A-Za-z]", x))))
    return df


def init_motifold_df(fasta_file: str, predictions_file: str) -> pd.DataFrame:
    loader = DataLoader(fasta_file, predictions_file)
    df = loader.build_dataframe()
    df = (
        df
        .pipe(parse_id_column)
        .pipe(reorder_columns)
        .pipe(adjust_mfe_to_kcal)
        .pipe(rename_columns)
        .pipe(fill_nan_motifs)
        .pipe(add_relative_start_end)
        .pipe(add_distance_to_mfe)
        .pipe(add_clip_window, CLIP_NT_WINDOW)
        .pipe(add_motifs_clip_window)
        .pipe(replace_ambiguous_motifs, AMBIGUOUS_MOTIFS)
    )
    return df

# TODO: Implement method.
def init_motices_df(fasta_file: str, predictions_file: str) -> pd.DataFrame:
    raise NotImplementedError("init_motices_df method is not implemented yet.")