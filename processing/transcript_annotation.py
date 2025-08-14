"""
transcript_annotation.py
~~~~~~~~~~~~~~~~~~~~~~~~

Includes all classes related to the annotation and parsing of the GFF3 files.

Classes
-------
TranscriptAnnotation
    Class representing an entry in the GFF3 file with its related metadata.
TranscriptAnnotationParser
    Parses, validates and processes the entries in the GFF3 file.

author: U.B.
"""

from dataclasses import dataclass, field
from intervaltree import IntervalTree

from processing.clip_annotation import ClipEntry
from processing.utils import normalize_chr

@dataclass(slots=True)
class TranscriptAnnotation:
    """
    TranscriptAnnotation class representing an entry in a GFF3 file.

    Attributes
    ----------
    chromosome: str
        The chromosome the sequence of the annotation is located in.
    start: str
        Start position of the sequence.
    end: str
        End position of the sequence.
    strand_orientation: str
        The orientation in the strand of the sequence.
    feature_type: str
        The feature type of the sequence.
    transcript_attributes: dict[str, str]
        The attributes of the annotation.
    clip_data : list[ClipEntry]
        Unique CLIP entries that overlap with the transcript position.
        Internally stored as a set to enforce uniqueness.
    sequence: str
        DNA sequence of the transcript.
    """
    chromosome: str
    start: int
    end: int
    strand_orientation: str
    feature_type: str
    transcript_attributes: dict[str, str]
    _clip_data: set[ClipEntry] = field(default_factory=set, init=False, repr=False)
    sequence: str = ""

    @property
    def clip_data(self) -> list[ClipEntry]:
        """
        Unique CLIP entries that overlap with the transcript position.
        Returned as a list for compatibility, but internally stored as a set
        to ensure uniqueness.

        Returns
        -------
        list[ClipEntry]
            The set of CLIP entries as a list.
        """
        return list(self._clip_data)

    def to_dict(self) -> dict[str, str | int| dict | list]:
        """
        Transforms the value of the attributes into a dictionary and returns it.

        Returns
        -------
        dict[str, str | int| dict | list]
            The transformed dictionary.
        """
        return {
            "chromosome": self.chromosome,
            "start": self.start,
            "end": self.end,
            "strand orientation": self.strand_orientation,
            "feature type": self.feature_type,
            "attributes": self.transcript_attributes,
            "sequence": self.sequence,
            "clip data": [clip.to_dict() for clip in self.clip_data]
        }
    
    def add_clip(self, clip: ClipEntry) -> None:
        """
        Adds a CLIP entry to the set of CLIP data.

        Parameters
        ----------
        clip: ClipEntry
            CLIP entry to add to the clip_data set.
        """
        self._clip_data.add(clip)

class TranscriptAnnotationParser:
    """
    TranscriptAnnotationParser class responsible for parsing the GFF3 file.

    Attributes
    ----------
    gff3_file: str
        The path to the GFF3 file.
    included_features: list[str]
        The feature types to be included when parsing.
    """
    def __init__(self, gff3_file: str, included_features: list[str]):
        """
        Initializes a TranscriptAnnotationParser instance.

        Parameters
        ----------
        gff3_file: str
            Path to the GFF3 file.
        included_features: list[str]
            The features to be included. 
            All other features are skipped.
        """
        self.gff3_file = gff3_file
        self.included_features = included_features
    
    def _validate_entry(self, line: str, chromosome: str) -> list[str] | None:
        """
        Validates the line and entry. Returns the valid entry as a list, otherwise returns none.

        Parameters
        ----------
        line: str
            Line in the GFF3 file with the entry information.
        chromosome: str
            Chromosome the transcripts are parsed from.

        Returns
        -------
        list[str]
            List with the data of the entry.
        """
        if line.startswith("#"):
            return None
        entry = line.strip().split("\t")
        if entry[2] not in self.included_features:
            return None
        if normalize_chr(entry[0]) != chromosome:
            return None
        return entry
    
    def _process_entry(self, entry: list[str], chromosome: str) -> tuple[int, int, str, TranscriptAnnotation] | None:
        """
        Processes the entry and returns a tuple of the start position, end position, transcript ID and the transcript annotation.
        If there is no transcript ID to get, it returns none.

        Parameters
        ----------
        entry: list[str]
            Entry of the GFF3 file as a list.
        chromosome: str
            Chromosome the transcript is located in.

        Returns
        -------
        tuple[int, int, str, TranscriptAnnotation] 
            Tuple containing:
            int: Starting position of the transcript.
            int: Ending position of the transcript.
            str: ID of the transcript.
            TranscriptAnnotation: Annotation of the transcript.
        """
        start, end = int(entry[3]) - 1, int(entry[4])  # GFF3 is 1-based
        attrs = self._parse_attributes(entry[8])
        transcript_id = attrs.get("ID")
        if not transcript_id:
            return None
        transcript_annotation = TranscriptAnnotation(
            chromosome=chromosome,
            start=start,
            end=end,
            strand_orientation=entry[6],
            feature_type=entry[2],
            transcript_attributes=attrs,
        )
        return start, end, transcript_id, transcript_annotation


    @staticmethod
    def _parse_attributes(attributes_str: str) -> dict[str, str]:
        """
        Parses the attributes of a GFF3 entry into a dictionary and returns it

        Parameters
        ----------
        attributes_str: str
            String containing all atributes. The attributes are separated by a ';'
            and the values are assigned with a '=' to their keys.

        Returns
        -------
        dict[str, str]
            All attributes as key-value pairs.
        """
        attributes = {}
        for attribute in attributes_str.strip().split(";"):
            if "=" in attribute:
                key, val = attribute.split("=", 1)
                attributes[key] = val
        return attributes
    
    def build_interval_tree_for_chrom(self, chromosome: str) -> tuple[dict[str, TranscriptAnnotation], IntervalTree]:
        """
        Builds an intervaltree for the specified chromosome, where the start and end of the corresponding transcript is stored
        and creates a dictionary with the transcripts.

        Parameters
        ----------
        chromosome: str
            Specified chromosome the transcripts are located in.

        Returns
        -------
        tuple[dict[str, TranscriptAnnotation], IntervalTree]
            Tuple containing:
            dict[str, TranscriptAnnotation]: Diciotnary of the
            transcripts with their IDs and the annotation.
            IntervalTree: The IntervalTree with the stored starting position, 
            ending position and the correpsondint transcript ID.
        """
        transcripts, tree = {}, IntervalTree()
        with open(self.gff3_file) as f:
            for line in f:
                entry = self._validate_entry(line, chromosome)
                if entry is None:
                    continue
                result: tuple[int, int, str, TranscriptAnnotation] = self._process_entry(entry, chromosome)
                if result is None:
                    continue
                start, end, transcript_id, transcript_annotation = result
                transcripts[transcript_id] = transcript_annotation
                tree.addi(start, end, transcript_id)
        return transcripts, tree

    @staticmethod
    def _is_valid_chromosome(chromosome: str) -> bool:
        """
        Check if chromosome name is shorter than six characters.
        
        Parameters
        ----------
        chromosome: str
            The chromosome to validate.

        Returns
        -------
        bool
            True if the chromosome has a length shorter than six characters.
            False if the chromosome is longer than five characters.
        """
        return len(chromosome) <= 6

    @staticmethod
    def _parse_chromosome(line: str) -> str | None:
        """
        Extract and normalize chromosome name from a GFF3 line,
        or None if the line is a comment.

        Parameters
        ----------
        line: str
            Line with the entry where the chromosome is parsed from.

        Returns
        -------
        str | None
            The normalized chromosome if the line is an entry,
            otherwise None.
        """
        if line.startswith("#"):
            return None
        entry = line.strip().split("\t")
        return normalize_chr(entry[0])

    def extract_chromosomes(self) -> list[str]:
        """
        Extracts all valid chromosomes in the GFF3 file.

        Returns
        -------
        list[str]
            List with all chromosomes.
        """
        chromosomes = set()
        with open(self.gff3_file) as file:
            for line in file:
                chromosome = self._parse_chromosome(line)
                if chromosome and self._is_valid_chromosome(chromosome):
                    chromosomes.add(chromosome)
        return list(chromosomes)