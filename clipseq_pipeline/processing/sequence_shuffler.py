import random

import numpy as np
from numpy.typing import NDArray
from dinuc_shuf import shuffle

def mono_nt_preserve_shuffler(sequence: str) -> str:
    """
    Shuffles the nucleotides of a sequence by preserving the mononucleotides,
    does not preserve dinucleotide.

    Parameters
    ----------
    sequence: str
        The string of the sequence to shuffle the nucleotides.

    Returns
    -------
    str
        The string of the shuffled sequence.
    """
    sequence_list = list(sequence)
    random.shuffle(sequence_list)
    return "".join(sequence_list)

def di_nt_preserve_shuffler(sequence: str) -> str:
    """
    Generates randomized sequences that preserve the dinucleotide composition
    of the original input using a one-hot-encoded dinucleotide shuffle.

    Parameters
    ----------
    sequence: str
        The string of the sequence to shuffle the nucleotides.

    Returns
    -------
    str
        The string of the shuffled sequence.
    """
    seq_alphabet = np.array(["A", "C", "G", "T"], dtype="S1")

    def one_hot_encode(sequence: str, dtype: type[np.uint8] = np.uint8) -> NDArray[np.uint8]:
        """Convert A/C/G/T string to one-hot array of shape (L, 4)."""
        seq = sequence.upper()
        seq_arr = np.frombuffer(seq.encode("utf8"), dtype="S1")
        return (seq_arr[:, None] == seq_alphabet[None, :]).astype(dtype)

    def one_hot_decode(one_hot: NDArray[np.uint8]) -> str:
        """Convert one-hot array back to A/C/G/T string."""
        return seq_alphabet[one_hot.argmax(axis=1)].tobytes().decode("utf8")

    one_hot = one_hot_encode(sequence)
    shuffled_oh = shuffle(one_hot[None, :, :])
    shuffled = one_hot_decode(shuffled_oh[0])
    return shuffled