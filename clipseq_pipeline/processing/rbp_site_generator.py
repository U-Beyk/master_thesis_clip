'''
rbp_site_generator.py
=====================

Module containing the class for generating RNA-binding-protein binding sites.

author: U.B.
'''

from collections.abc import Iterator
import os
from typing import Callable

from ..constants import RBP_NT_LENGTH
from .clip_processing import ClipEntry, ClipProcessor
from .genome_processor import GenomeProcessor
from .gff3_processor import Gff3Processor
from .sequence_shuffler import di_nt_preserve_shuffler

type SequenceFn = Callable[[ClipEntry], str]
    
class RbpSiteGenerator:
    '''
    Class generating RNA-binding-protein binding sites for each CLIP entry.

    Attributes
    ----------
    organism: str
        String of the organism name.
    clip_data: ClipProcessor
        ClipProcessor instance, where the CLIP data can be filtered from.
    genome: GenomeProcessor
        GenomeProcessor instance to filter the DNA sequence of the rbp sites.
    gff3_index: Gff3Processor
        Gff3Processor instance to get the feature types and annotation
        information of the rbp sites.
    '''

    def __init__(self, organism: str):
        '''
        Initializes a RbpSiteGenerator instance.

        Parameters
        ----------
        organism: str
            Name string of the organism, the rbp sites are from.
        '''
        self.organism = organism
        self.clip_data = ClipProcessor(f"data/datasets/{self.organism}/{self.organism}_clip.txt")
        self.genome = GenomeProcessor(f"data/datasets/{self.organism}/{self.organism}_genome.fa")
        self.gff3_index = Gff3Processor(f"data/datasets/{self.organism}/{self.organism}_annotations.gff3")

    def _write_fasta_file(self, filename: str, sequence_func: SequenceFn) -> None:
        """
        Writes CLIP entries to a FASTA file using the provided sequence function.

        Parameters
        ----------
        filename: str
            Name of the output FASTA file.
        sequence_func: SequenceFn
            Function that takes a clip and returns a string
            representation for the FASTA file.
        """
        output_dir = f"./data/fasta_files/{self.organism}"
        os.makedirs(output_dir, exist_ok=True)
        filepath = os.path.join(output_dir, filename)

        with open(filepath, "w") as file:
            for clip_counter, clip in enumerate(self._iterate_rbpsites(), start=1):
                clip.clip_id = f"{self.organism}:{clip_counter}"
                file.write(sequence_func(clip))

    def write_fasta(self) -> None:
        '''
        Assigns an ID to the CLIP entries and writes them into a FASTA file.
        '''
        self._write_fasta_file(
            f"{self.organism}_rbp_sites.fasta",
            lambda clip: clip.to_fasta()
        )

    def write_shuffled_fasta(self) -> None:
        '''
        Assigns an ID to the CLIP entries and writes them with their
        dinucleotide shuffled sequence into a FASTA file.
        '''
        self._write_fasta_file(
            f"{self.organism}_shuffled_rbp_sites.fasta",
            lambda clip: clip.to_shuffled_fasta()
        )

    def _iterate_rbpsites(self) -> Iterator[ClipEntry]:
        '''
        Iterates through the CLIPs and yields the entries,
        to which features can be assigned,
        to which a sequence can be assigned
        and that is unique.

        Yields
        ------
        ClipEntry
            CLIP with a sequence and features that consequently correspond to a rbp site.
        '''
        unique_clips = set()
        for clip in self.clip_data.iterate_clips():
            if not self._assign_features(clip):
                continue
            if not self._assign_sequence(clip):
                continue
            if clip in unique_clips:
                continue
            unique_clips.add(clip)
            yield clip

    def _assign_features(self, clip: ClipEntry) -> bool:
        '''
        Assigns the feature types and returns true if the CLIP entry has features.

        Parameters
        ----------
        clip: ClipEntry
            The CLIP entry to check features of and assign them.

        Returns
        -------
        bool
            True if the CLIP entry has features.
        '''
        features = []
        if self.gff3_index:
            features = self.gff3_index.get_features(
                clip.chromosome, clip.strand_orientation, clip.clip_start, clip.clip_end
            )

        if not features:
            return False
        else:
            clip.feature_types = features
            return True

    def _assign_sequence(self, clip: ClipEntry) -> bool:
        """
        Checks if a CLIP entry has a sequence and assigns it to the entry.
        Also assigns the sequence start and end to the CLIP entry.
        And creates a dinucelotide preserving shuffle and a mononucleotide
        shuffle of the sequence.

        Parameters
        ----------
        clip: ClipEntry
            The CLIP entry to check a sequence of and assign it.

        Returns
        -------
        bool
            True if the CLIP entry has a sequence.
        """
        if not self.genome:
            return False
        
        desired_length = RBP_NT_LENGTH
        seq_start, seq_end = self._compute_extended_coordinates(clip, desired_length)
        seq_start, seq_end = self._adjust_for_chromosome_bounds(clip.chromosome, seq_start, seq_end, desired_length)
        sequence = self._fetch_sequence(clip, seq_start, seq_end)
        shuffled_sequence = di_nt_preserve_shuffler(sequence)
        
        if (seq_end - seq_start) < desired_length:
            return False
        else:
            clip.sequence_start, clip.sequence_end, clip.sequence = seq_start, seq_end, sequence
            clip.shuffled_sequence = shuffled_sequence
            return bool(sequence)
    
    def _compute_extended_coordinates(self, clip: ClipEntry, desired_length: int) -> tuple[int, int]:
        '''
        Compute symmetric extension around the clip to reach desired length.
        Considers the CLIP being near the start of the chromosome.

        Parameters
        ----------
        clip: ClipEntry
            ClipEntry to compute the extended coordinates of.
        desired_length: int
            The length the RNA sequence will be extended to.

        Returns
        -------
        tuple[int, int]
            Tuple consisting of the sequence start and sequence end.
        '''
        clip_start, clip_end = clip.clip_start, clip.clip_end
        clip_length = clip_end - clip_start

        extend = (desired_length - clip_length) // 2
        extend_left = min(extend, clip_start)
        sequence_start = clip_start - extend_left
        extend_right = desired_length - clip_length - extend_left
        sequence_end = clip_end + extend_right

        return sequence_start, sequence_end

    def _adjust_for_chromosome_bounds(
            self, 
            chromosome: str, 
            sequence_start: int, 
            sequence_end: int, 
            desired_length: int
    ) -> tuple[int, int]:
        '''
        Ensure extended region stays within chromosome bounds.

        Parameters
        ----------
        chromosome: str
            String of the chromosome name.
        sequence_start: int
            Start of the sequence.
        sequence_end: int
            End of the sequence.
        desired_length: int
            Length of the RNA sequence to adjust the start and end to.
        '''
        chr_len = self.genome.get_chromosome_length(chromosome)

        if sequence_end > chr_len:
            excess = sequence_end - chr_len
            shift_left = min(excess, sequence_start)
            sequence_start -= shift_left
            sequence_end = min(sequence_start + desired_length, chr_len)

        sequence_start = max(0, sequence_start)
        sequence_end = max(sequence_start, sequence_end)
        return sequence_start, sequence_end

    def _fetch_sequence(self, clip: ClipEntry, sequence_start: int, sequence_end: int) -> str:
        '''
        Fetch sequence from genome.

        Parameters
        ----------
        clip: ClipEntry
            ClipEntry to fetch the sequence of.
        sequence_start: int
            Start of the sequence (not the same as start of CLIP).
        sequence_end: int
            End of the sequence (not the same as end of CLIP).

        Returns
        -------
        str
            String of the sequence.
        '''
        return self.genome.get_sequence(
            clip.chromosome,
            sequence_start,
            sequence_end,
            clip.strand_orientation
        ) or ""