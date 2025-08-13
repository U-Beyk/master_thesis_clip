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

def normalize_chr(chromosome: str) -> str:
    chromosome = chromosome.upper()
    if chromosome.startswith("CHR"):
        chromosome = chromosome[3:]
    return chromosome