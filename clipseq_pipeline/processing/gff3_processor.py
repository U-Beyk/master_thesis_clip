'''
gff3_processor.py
=================

Module for parsing and processing GFF3 genomic annotation files.

author: U.B.
'''

from collections import defaultdict

from intervaltree import Interval, IntervalTree

from .utils import normalize_chr

class Gff3Processor:
    '''
    Class handling the processing and mapping of the GFF3 entries.

    Attributes
    ----------
    _interval_trees: dict[tuple[str, str], IntervalTree]
        Interval trees by chromosome and strand orientation.
    '''

    def __init__(self, gff3_file: str):
        '''
        Initializes a Gff3Processor object.

        Parameters
        ----------
        gff3_file: str
            String of the path to the gff3-file.
        '''
        self._interval_trees: dict[tuple[str, str], IntervalTree] = self._load_gff3(gff3_file)
    
    def _load_gff3(self, gff3_file: str) -> dict[tuple[str, str], IntervalTree]:
        '''
        Loads the gff3 file and creates interval trees.

        Parameters
        ----------
        gff3_file: str
            The path to the gff3 file.

        Returns
        -------
        dict[tuple[str, str], IntervalTree]
            Dictionary with the interval trees of each chromosome, strand orientation tuple.
        '''
        interval_trees: dict[tuple[str, str], IntervalTree] = defaultdict(IntervalTree)
        with open(gff3_file) as file:
            for line in file:
                if line.startswith("#"):
                    continue
                chromosome, feature_type, start, end, strand_orientation = self._process_line(line)
                if strand_orientation not in ("+", "-"):
                    continue
                interval_trees[(chromosome, strand_orientation)].addi(start - 1, end, feature_type)
        return interval_trees
    
    @staticmethod
    def _process_line(line: str) -> tuple[str, str, int, int, str]:
        '''
        Processes a line in the gff3 file and returns the necessary information to build an interval tree.

        Parameters
        ----------
        line: str
            A line from the gff3 file.

        Returns
        -------
        tuple[str, str, int, int, str]
            Tuple containing the chromosmoe string, the feature_type, the start, end and strand orientation
            of the annotation/line.
        '''
        columns = line.strip().split("\t")
        chromosome, _, feature_type, start, end, _, strand_orientation, _, _ = columns
        chromosome = normalize_chr(chromosome)
        start, end = int(start), int(end)
        return chromosome, feature_type, start, end, strand_orientation

    def get_features(self, chromosome: str, strand_orientation: str, start: int, end: int) -> list[str]:
        '''
        Gets the feature type of the gff3 annotation by checking for overlaps with the interval trees of a specific sequence.

        Parameters
        ----------
        chromosome: str
            The string of the chromosome.
        strand_orientation: str
            The strand orientation of the sequence.
        start: int
            Start of the sequence.
        end: int
            End of the sequence.

        Returns
        -------
        list[str]
            A list of feature type names that overlap the specified genomic interval.
            Returns an empty list if no overlapping features are found.
        '''
        interval_tree = self._interval_trees.get((chromosome, strand_orientation))
        if not interval_tree:
            return []
        hits: set[Interval] = interval_tree.overlap(start, end)
        return list({interval.data for interval in hits})