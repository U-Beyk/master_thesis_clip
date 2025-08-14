"""
genome_processing.py
~~~~~~~~~~~~~~~~~~~~

This file contains the class for handling and processing the Genome of a specified organism.

Classes
-------
Genome
    Class handling the genome data.

author: U.B.
"""

from pyfaidx import Fasta, Sequence
from processing.utils import normalize_chr

class Genome:
    """
    Class for handling genome related data.

    Attributes
    ----------
    genome: Fasta
        The genome of the specified organism as a Fasta object.
    _chromosome_name_map: dict[str, str]
        Dictionary mapping the chromosome to its normalized variant. 
        This attribute is private.
    """
    def __init__(self, genome_fasta: str):
        """
        Initializes a Genome object corresponding to the genome of the specified organism.

        Parameters
        ----------
        genome_fasta: str
            The path to the fasta file of the specified organism.
        """
        self.genome = Fasta(genome_fasta)
        self._chromosome_name_map: dict[str, str] = {
            normalize_chr(chromosome): chromosome 
            for chromosome in self.genome.keys()
            }

    def _map_chromosome_name(self, chromosome: str) -> str | None:
        """
        Return the actual chromosome name from the FASTA, or None if not found.
        
        Parameters
        ----------
        chromosome (str): 
            The specified chromosome name/number as a string.

        Returns
        -------
        str | None
            The mapped chromosome, if there are any, otherwise None.
        """
        return self._chromosome_name_map.get(normalize_chr(chromosome))

    def _get_raw_sequence(self, actual_chromosome: str, start: int, end: int) -> str:
        """
        Fetch raw sequence from the genome.

        Parameters
        ----------
        actual_chromosome: str
            The chromosome the sequence have to be fetched from.
        start: int
            Starting base of the sequence.
        end: int
            Ending base of the sequence.

        Returns
        -------
        str
            The raw sequence.
        """
        sequence_record: Sequence = self.genome[actual_chromosome][start:end]
        return sequence_record.seq
    
    @staticmethod
    def _reverse_complement(sequence: str) -> str:
        """
        Return the reverse complement of a DNA sequence.

        Parameters
        ----------
        sequence: str
            The DNA sequence.

        Returns
        -------
        str
            The reverse complement of the sequence.
        """
        complement_table = str.maketrans("ACGTacgt", "TGCAtgca")
        return sequence.translate(complement_table)[::-1]

    def extract_sequence(self, chromosome: str, start: int, end: int, strand_orientation: str) -> str:
        """
        Get strand-aware genomic sequence, without ambiguous nucleotides.

        Parameters
        ----------
        chromosome: str
            The chromosome to extract the sequence from.
        start: int
            Starting base of the sequence.
        end: int
            Ending base of the sequence.
        strand_orientation: str
            The strand orientation of the sequence.

        Returns
        -------
        str
            The DNA sequence in the genome.
        """
        actual_chrom = self._map_chromosome_name(chromosome)
        if not actual_chrom:
            return ""
        sequence = self._get_raw_sequence(actual_chrom, start, end)
        if any(nucleotide not in set("ACGTacgt") for nucleotide in sequence):
            return ""
        return self._reverse_complement(sequence) if strand_orientation == "-" else sequence