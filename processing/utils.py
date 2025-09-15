"""
utils.py
~~~~~~~~

This module contains generic utility functions.

Functions
---------
normalize_chr
    Normalizes the chromosome string.

author: U.B.
"""

# Maps ambiguous chromosome names
AMBIGUOUS_CHR_MAPPING = {
    # Mitochondrial variants
    "M": "MT",
    "MTDNA": "MT",
    "MITO": "MT"
}

def normalize_chr(chromosome: str) -> str:
    '''
    Normalizes the chromosome, by turning it into uppercase, cutting of the "CHR" string at the beginning 
    if it exists and standardizes ambiguous chromosome names.

    Parameters
    ----------
    chromosome: str
        The chromosome name as a string.

    Returns
    -------
    str
        The normalized chromosome name.
    '''
    chromosome = chromosome.upper()
    if chromosome.startswith("CHR"):
        chromosome = chromosome[3:]
    return AMBIGUOUS_CHR_MAPPING.get(chromosome, chromosome)