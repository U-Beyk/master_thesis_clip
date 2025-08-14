"""
clip_annotation.py
~~~~~~~~~~~~~~~~~~

This file contains all CLIP-related classes for processing and annotating CLIP data.

Classes
-------
ClipEntry
    Represents a blueprint for clip data entries with attributes like binding site coordinates,
    strand orientation, RNA binding protein name, and related metadata.
ClipParser
    Parses the .txt file (by chromosome) containing the CLIP data

author: U.B.
"""

from dataclasses import dataclass
from processing.utils import normalize_chr

@dataclass(frozen=True, slots=True)
class ClipEntry:
    """
    ClipEntry class representing the blueprint for clip data entries.

    Attributes
    ----------
    start: int 
        Start of the binding site.
    end: int
        End of the binding site.
    strand_orientation: str
        Strand orientation of the sequence.
    rbp_name: str
        Name of the RNA binding protein.
    methods: str
        Experiment method.
    software: str
        Peak calling software.
    sample: str
        Sample/tisue used in study.
    accession_data: str
        Accession of raw data.
    confidence_score: float
        Confidence score assigned by computational method.
    """
    start: int
    end: int
    strand_orientation: str
    rbp_name: str
    method: str
    software: str
    sample: str
    accession_data: str
    confidence_score: float

    def to_dict(self) -> dict[str, str | int| float]:
        """
        Transforms the attributes into a dictionary and returns it.

        Returns
        -------
        dict[str, str | int | float]
            Dictionary representation of the attributes.
        """
        return {
            "start": self.start,
            "end": self.end,
            "strand orientation": self.strand_orientation,
            "RNA binding protein": self.rbp_name,
            "method": self.method,
            "software": self.software,
            "sample/tissue": self.sample,
            "accession of raw data": self.accession_data,
            "confidence score": self.confidence_score
        }
    
class ClipParser:
    """
    ClipParser class responsible for the parsing of the .txt file containing the CLIP data.

    Attributes
    ----------
    clip_file: str
        Path of the CLIP file.
    """
    def __init__(self, clip_file: str):
        """
        Initializes a ClipParser instance.

        Parameters
        ----------
        clip_file: str
            Path to the .txt file with the CLIP data.
        """
        self.clip_file = clip_file

    def filter_clips_of_chromosome(self, chromosome: str) -> list[ClipEntry]:
        """
        Filters the CLIP data by the specified chromosome and returns it as a list.

        Parameters
        ----------
        chromosome: str
            Chromosome number/name.

        Returns
        -------
        list[ClipEntry]
            All CLIP entries in the specified chromosome.
        """
        clip_data = []
        with open(self.clip_file) as f:
            for line in f:
                entry = line.strip().split("\t")
                clip_chrom = normalize_chr(entry[0])
                if clip_chrom != chromosome:
                    continue
                clip_data.append(self._parse_clip_entry(entry))
        return clip_data

    @staticmethod
    def _parse_clip_entry(entry: list[str]) -> ClipEntry:
        """
        Convert a line's fields into a ClipEntry object with correct types.
        
        Parameters
        ----------
        entry: list[str]
            CLIP entry.

        Returns
        -------
        ClipEntry
            ClipEntry object of the line corresponding to an entry.
        """
        # The column with the method and software sometimes only contain one entry. 
        # In those cases the method and software attributes are set the same.
        method_software = entry[6].split(",")
        method, software = (method_software if len(method_software) == 2 
                            else (method_software[0], method_software[0]))
        return ClipEntry(
            start=int(entry[1]),
            end=int(entry[2]),
            strand_orientation=entry[4],
            rbp_name=entry[5],
            method=method,
            software=software,
            sample=entry[7],
            accession_data=entry[8],
            confidence_score=float(entry[9])
        )