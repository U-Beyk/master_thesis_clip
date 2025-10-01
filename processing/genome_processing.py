"""
genome_processing.py
~~~~~~~~~~~~~~~~~~~~

This file contains the class for handling and processing the Genome.

Classes
-------
GenomeProcessing
    Class handling the genome data.

author: U.B.
"""

from pyfaidx import Fasta
from processing.utils import normalize_chr

class GenomeProcessing:
    '''
    Class handling and processing genome data.

    Attributes
    ----------
    genome: Fasta
        The genome as a Fasta object.
    _chromosome_name_map: dict[str, str]
        Mapping of the normalized chromosome to the chromosome names in the genome fasta.
    '''

    def __init__(self, genome_fasta: str):
        '''
        Initializes a GenomeProcessing object to extract sequence data.

        Parameters
        ----------
        genome_fasta: str
            Path to the fasta file of the genome.
        '''
        self.genome = Fasta(genome_fasta)
        self._chromosome_name_map: dict[str, str] = self._map_chromosome_to_name()
    
    def _map_chromosome_to_name(self) -> dict[str, str]:
        '''
        Maps the normalized chromosomes to the chromosome names in the fasta file of the genome.

        Returns
        -------
        dict[str, str]
            The dictionary with the mapped chromosomes.
        '''
        chromosome_name_map = dict()
        for chromosome in self.genome.keys():
            if len(chromosome) <= 6:
                chromosome_name_map[normalize_chr(chromosome)] = chromosome
        return chromosome_name_map

    def get_sequence(self, chromosome: str, start: int, end: int, strand_orientation: str) -> str | None:
        '''
        Gets the sequence from the genome by specifiying the chromosome, start, end and strand orientation of the sequence.

        Parameters
        ----------
        chromosome: str
            The chromosome the sequence is located in.
        start: int
            Start of the sequence.
        end: int
            End of the sequence.
        strand_orientation: str
            The orientation of the strand. Either "+" or "-", anything else will be ignored.

        Returns
        -------
        str | None:
            The specified sequence or None, if the sequence contains ambiguous nucleotides.
        '''
        chromosome = self._chromosome_name_map.get(normalize_chr(chromosome))
        sequence = str(self.genome[chromosome][start:end]).upper()
        if not set(sequence).issubset({"A", "C", "G", "T"}):
            return None
        if strand_orientation == "-":
            return self.reverse_complement(sequence)
        elif strand_orientation == "+":
            return sequence
        else:
            return None

    @staticmethod
    def reverse_complement(sequence: str) -> str:
        '''
        Forms the reverse complement of the sequence.

        Parameters
        ----------
        sequence: str
            The sequence string to form a reverse complement of.
        
        Returns
        -------
        str
            The reverse complement of the sequence.
        '''
        complement = str.maketrans("ACGTacgt", "TGCAtgca")
        return sequence.translate(complement)[::-1]