'''
rbp_site_generator.py
~~~~~~~~~~~~~~~~~~~~~

This file contains the class for generating RNA binding protein binding sites.

Classes
-------
RbpSiteGenerator
    Class for generating rbp sites.

author: U.B.
'''

from collections.abc import Iterator

from processing.clip_processing import ClipEntry, ClipProcessing
from processing.genome_processing import GenomeProcessing
from processing.gff3_processing import Gff3Processing
    
class RbpSiteGenerator:
    '''
    Class generating RNA binding protein binding sites for each CLIP entry.

    Attributes
    ----------
    organism: str
        String of the organism name.
    clip_data: ClipProcessing
        ClipProcessing instance, where the CLIP data can be extracted.
    genome: GenomeProcessing
        GenomeProcessing instance to filter the DNA sequence of the rbp sites.
    gff3_index: Gff3Processing
        Gff3Processing inctance to get the feature types and annotation
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
        self.clip_data = ClipProcessing(f"data/datasets/{self.organism}/{self.organism}_clip.txt")
        self.genome = GenomeProcessing(f"data/datasets/{self.organism}/{self.organism}_genome.fa")
        self.gff3_index = Gff3Processing(f"data/datasets/{self.organism}/{self.organism}_annotations.gff3")

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
        clip.feature_types = features
        return True

    def _assign_sequence(self, clip: ClipEntry) -> bool:
        '''
        Checks if a CLIP entry has a sequence and assigns it to the entry.
        Also assigns the sequence start and end to the CLIP entry.

        Parameters
        ----------
        clip: ClipEntry
            The CLIP entry to check a sequence of and assign it.

        Returns
        -------
        bool
            True if the CLIP entry has a sequence.
        '''
        if not self.genome:
            return False
        seq_start = max(1, clip.clip_start - 50)
        seq_end = clip.clip_end + 50
        clip.sequence_start, clip.sequence_end  = seq_start, seq_end
        clip.sequence = self.genome.get_sequence(
            clip.chromosome, seq_start, seq_end, clip.strand_orientation
        ) or ""
        return bool(clip.sequence)

    def iterate_rbpsites(self) -> Iterator[ClipEntry]:
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

    def write_fasta(self) -> None:
        '''
        Assigns an ID to the CLIP entries and writes them into a FASTA file.
        '''
        with open(f"data/fasta_files/{self.organism}_rbp_sites.fasta", "w") as file:
            for clip_counter, clip in enumerate(self.iterate_rbpsites(), start=1):
                clip.clip_id = f"{self.organism}:{clip_counter}"
                clip_string = clip.to_fasta()
                file.write(clip_string)