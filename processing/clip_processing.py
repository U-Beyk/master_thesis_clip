"""
clip_annotation.py
~~~~~~~~~~~~~~~~~~

This file contains all CLIP-related classes for processing and annotating CLIP data.

Classes
-------
ClipEntry
    Represents a blueprint for clip data entries with attributes like binding site coordinates,
    strand orientation, RNA binding protein name, and related metadata.
ClipProcessing
    Parses the .txt file and processes the CLIP data.

author: U.B.
"""

from collections import defaultdict
from collections.abc import Iterator
from dataclasses import dataclass

from processing.utils import normalize_chr

# List of peak-calling tools where higher scores indicate higher-confidence peaks.
# MOVE IT TO THE CONFIG??????????
HIGH_CONFIDENCE_SCORER = [
    "eCLIP", "PARalyzer", "PureCLIP", "CTK", "CIMS", "MiClip",
    "Piranha_0.01", "PIP-seq"
]

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
        Sample/tissue used in the expeiment.
    accession_data: str
        Raw accession data.
    accesson_experiment: str
        String of the experiment.
    confidence_score: float
        Score of the CLIP entry.
    feature_types: list[str]
        List of feature types the CLIP entry corresponds to
        in the genome.
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
    sequence_start: int
    sequence_end: int

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
            self.method,
            self.software,
            self.sample,
            self.accession_data,
            self.accession_experiment,
            self.confidence_score
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
            self.method == other.method and
            self.software == other.software and
            self.sample == other.sample and
            self.accession_data == other.accession_data and
            self.accession_experiment == other.accession_experiment and
            self.confidence_score == other.confidence_score
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
    
class ClipProcessing:
    '''
    Class handling the parsing and processing of the Clip data.

    Attributes
    ----------
    clip_file: str
        Path to the clip file.
    _number_software_scores_experiments: dict[str, dict[str, float]]
        Scores of the tools in unique experiments corresponding to the top 10%,
        used to cut off data, that is below these scores. 
    '''

    def __init__(self, clip_file: str):
        '''
        Initializes a ClipProcessing object.

        Parameters
        ----------
        clip_file: str
            The path to the CLIP file.
        '''
        self.clip_file = clip_file
        self._number_software_scores_experiments= self._filter_unique_experiments()

    def _read_file(self) -> Iterator[list[str]]:
        '''
        Read the CLIP file line by line and yield parsed entries.

        Yields
        ------
        list[str]
            A list of fields from each line in the file, split by tab.
        '''
        with open(self.clip_file) as file:
            for line in file:
                yield line.strip().split("\t")

    def _parse_entry_for_filtering(self, entry: list[str]) -> tuple[str, str, float]:
        '''
        Parses an entry for the filtering step, where only the experiment, software and score are important.

        Parameters
        ----------
        entry: list[str]
            The content of an entry as a list of strings

        Returns
        -------
        tuple[str, str, float]
            Tuple containing the experiment string, software string and the score number as a float.
        '''
        score = float(entry[9])
        accession_parts = entry[8].split(",")
        experiment = accession_parts[1] if len(accession_parts) == 2 else accession_parts[0]
        method_parts = entry[6].split(",")
        software = method_parts[1] if len(method_parts) == 2 else method_parts[0]
        return experiment, software, score

    def _collect_scores(self) -> dict[str, dict[str, list[float]]]:
        '''
        Creates a dictioanry of scores for each experiment and the software tools used in that expeiment.

        Returns
        -------
        dict[str, dict[str, list[float]]]
            Dictionary of experiments and software tools with their scores
        '''
        experiments = defaultdict(lambda: defaultdict(list))
        for entry in self._read_file():
            experiment, software, score = self._parse_entry_for_filtering(entry)
            experiments[experiment][software].append(score)
        return experiments

    def _calculate_cutoff(self, scores: list[float], software: str) -> float:
        '''
        Calculates the cutoff to filter the CLIP data of.

        Parameters
        ----------
        scores: list[float]
            List of scores of the specified software tool.
        software: str
            Used software tool to calculate the cutoff of.

        Retruns
        -------
        float
            The cutoff score.
        '''
        scores.sort(reverse=True)
        decile_index = max(1, int(len(scores) * 0.1)) - 1  # ensure valid index
        if software in HIGH_CONFIDENCE_SCORER:
            return scores[decile_index]
        else:
            return scores[-decile_index]

    def _filter_unique_experiments(self) -> dict[str, dict[str, float]]:
        '''
        Filters the unique experiments with their software tools and the corresponding scores to cut off the CLIP data.

        Returns
        -------
        dict[str, dict[str, float]]
            Dictionary of experiments and software tools with their score cut offs.
        '''
        collected_scores = self._collect_scores()
        final_experiments = {}
        for experiment, software_dict in collected_scores.items():
            final_experiments[experiment] = {
                software: self._calculate_cutoff(scores, software)
                for software, scores in software_dict.items()
            }
        return final_experiments

    def iterate_clips(self) -> Iterator[ClipEntry]:
        '''
        Read the CLIP file line by line and yield the processed CLIP entry.
        Also filters the data by cutting off CLIPs that are too long 
        and have a too low score.

        Yields
        ------
        ClipEntry
            Entry of the CLIP data.
        '''
        for entry in self._read_file():
            # Parses the entry
            chromosome, start, end = normalize_chr(entry[0]), int(entry[1]), int(entry[2])
            strand, rbp_name = entry[4], entry[5]
            method_parts = entry[6].split(",")
            method, software = method_parts[0], method_parts[-1]
            sample = entry[7]
            acc_parts = entry[8].split(",")
            accession_data, accession_experiment = acc_parts[0], acc_parts[-1]
            score = float(entry[9])
            cutoff = self._number_software_scores_experiments.get(accession_experiment, {}).get(software)
            # Checks if the entry is valid
            if (end - start) > 100:
                continue
            if cutoff is None:
                continue
            if (software in HIGH_CONFIDENCE_SCORER and score < cutoff) or \
               (software not in HIGH_CONFIDENCE_SCORER and score > cutoff):
                continue
            # Creates and yields a ClipEntry object
            yield ClipEntry(
                clip_id="", chromosome=chromosome, clip_start=start, clip_end=end, 
                strand_orientation=strand, rbp_name=rbp_name, method=method, software=software, 
                sample=sample, accession_data=accession_data, accession_experiment=accession_experiment,
                confidence_score=score, feature_types=[], sequence="", sequence_start=None, sequence_end=None
            )