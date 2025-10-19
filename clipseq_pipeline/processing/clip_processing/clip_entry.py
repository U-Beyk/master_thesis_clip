'''
clip_entry.py
=============

This module contains the ClipEntry class, which corresponds to a single entry of the CLIP data used.

author: U.B.
'''

from dataclasses import dataclass

@dataclass(slots=True)
class ClipEntry:
    '''
    Class representing a CLIP entry.

    Attributes
    ----------
    clip_id: str
        unique ID string for the CLIP entry.
    chromosome: str
        String of the chromosome name/number
    clip_start: int
        Start of the CLIP entry.
    clip_end: int
        End of the CLIP entry.
    strand_orientation: str
        Strand orientation of the CLIP.
    rbp_name: str
        Name of the RNA binding protein.
    method: str
        CLIP method used in the experiment.
    software: str
        Software used for peak calling.
    sample: str
        Sample/tissue used in the experiment.
    accession_data: str
        Raw accession data.
    accesson_experiment: str
        String identifier of the experiment.
    confidence_score: float
        Score of the CLIP entry
        of the corresponding peak calling software.
    feature_types: list[str]
        List of feature types the CLIP entry corresponds to
        in the genome, like (gene, mRNA, lnc_RNA, etc.).
    sequence: str
        String of the DNA sequence.
    sequence_start: int | None
        Start of the sequence, if defined.
    sequence_end: int | None
        End of the sequence, if defined.
    '''

    clip_id: str
    chromosome: str
    clip_start: int
    clip_end: int
    strand_orientation: str
    rbp_name: str
    method: str
    software: str
    sample: str
    accession_data: str
    accession_experiment: str
    confidence_score: float
    feature_types: list[str]
    sequence: str
    sequence_start: int | None
    sequence_end: int | None

    def __hash__(self) -> int:
        '''
        Compute a hash value for a ClipEntry instance.

        Returns
        -------
        int
            The hash value representing this ClipEntry instance.
        '''
        return hash((
            self.chromosome,
            self.clip_start,
            self.clip_end,
            self.strand_orientation,
            self.rbp_name,
            self.sequence
        ))

    def __eq__(self, other: object) -> bool:
        '''
        Determine whether two ClipEntry instances are equal.

        Parameters
        ----------
        other: object
            The object to compare with this instance.

        Returns
        -------
        bool
            True if "other" is a ClipEntry instance and all attributes match,
            false otherwise.
        '''
        if not isinstance(other, ClipEntry):
            return NotImplemented
        return (
            self.chromosome == other.chromosome and
            self.clip_start == other.clip_start and
            self.clip_end == other.clip_end and
            self.strand_orientation == other.strand_orientation and
            self.rbp_name == other.rbp_name and
            self.sequence == other.sequence
        )

    def to_fasta(self) -> str:
        '''
        Creates the header and sequence of a CLIP entry.

        Returns
        -------
        str
            Header and sequence string.
        '''
        header = (
            f">id:{self.clip_id}|clip_range:{self.clip_start}-{self.clip_end}|rbp_name:{self.rbp_name}"
            f"|seq_range:{self.sequence_start}-{self.sequence_end}|features:{','.join(self.feature_types)}"
            )
        return f"{header}\n{self.sequence}\n"