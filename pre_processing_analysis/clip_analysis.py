"""
clip_analysis.py
~~~~~~~~~~~~~~~~

This module contains all classes for analyzing the CLIP data before processing.

Functions
---------
count_all_sequences_fasta
    Counts all sequences in all fasta files of a directory.
count_all_entries_gff3
    Counts all entries in a GFF3 file

Classes
-------
FullClipEntry
    Represents a blueprint for clip data entries, this includes the chromosome, peak_id
    and also the accession_experiment.
ClipDataAnalyzer
    Class responsible for the analysis of the CLIP data of a specified organism.
ClipSpeciesAnalyzer
    Analyzes all organisms specified in the cnfig file.

author: U.B.
"""

import matplotlib.pyplot as plt
import yaml
import os
from dataclasses import dataclass
from collections import Counter
from Bio import SeqIO

from processing.clip_processing import ClipEntry

# List of peak-calling tools where higher scores indicate higher-confidence peaks.
HIGH_CONFIDENCE_SCORER = [
    "eCLIP", "PARalyzer", "PureCLIP", "CTK", "CIMS", "MiClip",
    "Piranha_0.01", "PIP-seq"
]

def count_all_sequences_fasta(folder_path: str) -> int:
    """
    Counts all sequences in all fasta files in the specified directory.
    Prints the number of sequences for each fasta file.

    Parameters
    ----------
    folder_path: str
        String of the folder path where the fasta files are in.

    Returns
    -------
    int
        Total number of all sequences
    """
    total_sequences = 0
    for filename in os.listdir(folder_path):
        if filename.endswith(".fasta") or filename.endswith(".fa"):
            file_path = os.path.join(folder_path, filename)
            count = sum(1 for _ in SeqIO.parse(file_path, "fasta"))
            print(f"{filename}: {count} sequences")
            total_sequences += count
    return total_sequences

def count_all_entries_gff3(file_path: str) -> int:
    """
    Counts all entries in the specified GFF3 file.

    Parameters
    ----------
    file_path: str
        String of the file path.

    Returns
    -------
    int
        Total number of all entries.
    """
    unique_chromosomes = set()
    count = 0
    with open(file_path) as file:
        for line in file:
            if not line.startswith("#") and line.strip():
                entry = line.strip().split("\t")
                if len(entry[0]) < 7:
                    unique_chromosomes.add(entry[0])
                count += 1
    return count, unique_chromosomes

@dataclass(slots=True)
class FullClipEntry(ClipEntry):
    """
    FullClipEntry class representing the blueprint for clip data entries.
    This class includes atributes like chromosome, peak_id and accession_experiment.
    The rest is inherited from the ClipEntry class.

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
    peak_id: str
        ID of the peak in te POSTAR3 data.
    """
    chromosome: str
    peak_id: str

    def to_dict(self) -> dict[str, str | int| float]:
        """
        Transforms the attributes into a dictionary and returns it.

        Returns
        -------
        dict[str, str | int | float]
            Dictionary representation of the attributes.
        """
        return {
            "chromosome": self.chromosome,
            "start": self.clip_start,
            "end": self.clip_end,
            "peak ID": self.peak_id,
            "strand orientation": self.strand_orientation,
            "RNA binding protein": self.rbp_name,
            "method": self.method,
            "software": self.software,
            "sample/tissue": self.sample,
            "accession of raw data": self.accession_data,
            "accession experiment": self.accession_experiment,
            "confidence score": self.confidence_score
        }
    
class ClipDataAnalyzer:
    """
    Class for analysing the CLIP data.

    Attributes
    ----------
    organism: str
        Organism the CLIP data stems from.
    sowafter_tools: dict[str, bool]
        List of software tools to include in the analysis.
    clip_data: list[FullClipEntry]
        List of CLIP data.
    """

    __slots__ = ("organism", "clip_data", "software_tools")
    def __init__(self, organism: str, software_tools: list[str]):
        """
        Initializes a ClipDataAnalyzer instance.

        Parameters
        ----------
        organism: str
            Organism the CLIP data to extract from.
        software_tools: list[str]
            List of software tools to include in the analysis.
        """
        self.organism = organism
        self.software_tools = software_tools
        self.clip_data = self._extract_clip_data(organism)

    @property
    def unique_chromosomes(self) -> set[str]:
        """
        Gets and returns all unique chromosomes in the CLIP data.

        Returns
        -------
        set[str]
            Set of unique chromosome srings.
        """
        chromosomes = set()
        for data in self.clip_data:
            chromosomes.add(data.chromosome)
        return chromosomes

    @property
    def unique_methods(self) -> set[str]:
        """
        Gets and returns all unique methods in the CLIP data.

        Returns
        -------
        set[str]
            Set of unique methods strings.
        """
        methods = set()
        for data in self.clip_data:
            methods.add(data.method)
        return methods
    
    @property
    def unique_software_tools(self) -> set[str]:
        """
        Gets and returns all unique software tools.

        Returns
        -------
        set[str]
            Set of unique software tools strings.
        """
        tools = set()
        for data in self.clip_data:
            tools.add(data.software)
        return tools
    
    @property
    def unique_rbp(self) -> set[str]:
        """
        Gets and returns all unique RBPs.

        Returns
        -------
        set[str]
            Set of unique RBP strings.
        """
        tools = set()
        for data in self.clip_data:
            tools.add(data.rbp_name.upper())
        return tools
    
    @property
    def unique_accession_data(self) -> set[str]:
        """
        Gets and returns all unique accession entries.

        Returns
        -------
        set[str]
            Set of unique accession entries.
        """
        accession_data = set()
        for data in self.clip_data:
            accession_data.add(data.accession_data)
        return accession_data
    
    @property
    def unique_accession_experiments(self) -> set[str]:
        """
        Gets and returns all unique accession experiments.

        Returns
        -------
        set[str]
            Set of unique accession experiments.
        """
        accession_experiments = set()
        for data in self.clip_data:
            accession_experiments.add(data.accession_experiment)
        return accession_experiments
    
    def save_histogram(self, data: list[float], file_name: str) -> None:
        """
        Creates and stores a histogram for a specific list of data.

        Parameters
        ----------
        data: list[float]
            List of decimals/integers to be plotted.
        file_name: str
            String of th file name.
        """
        os.makedirs(f"plots/pre_processing_plots/{self.organism}", exist_ok=True)
        plt.figure(figsize=(8, 6))
        plt.hist(data, bins=50, edgecolor='black')
        plt.title(f"{file_name}")
        plt.xlabel("Value")
        plt.ylabel("Frequency")
        plt.savefig(f"plots/pre_processing_plots/{self.organism}/{file_name}.png", dpi=80)
        plt.close()

    def filter_entries(self, **criteria) -> list[FullClipEntry]:
        """
        Filters the data by the specified criterias.

        Parameters
        ----------
        **criteria: keyword arguments
            Arguments consisting of a key, value pair, the CLIP data is filtered by.
            The key corresponds to the attribute name and the value
            corresponds to the value of the attribute to filter by.

        Returns
        -------
        list[FullClipEntry]
            List of FullClipEntry data.
        """
        return [
            entry for entry in self.clip_data 
            if all(getattr(entry, attribute, None) == value for attribute, value in criteria.items())
            ]
    
    def all_data_scores_histogram(self) -> None:
        """
        Creates histograms for each unique software tool of their score distribution.
        """
        for tool in self.unique_software_tools:
            score_data = [entry.confidence_score for entry in self.filter_entries(software=tool)]
            if score_data:
                self.save_histogram(score_data, tool)

    def all_data_clip_length_histogram(self) -> None:
        '''
        Creates histograms for each unique method of the nucleotide length of the CLIP.
        '''
        for method in self.unique_methods:
            lengths = [entry.end - entry.start for entry in self.filter_entries(method=method)]
            if lengths:
                self.save_histogram(lengths, f"method_{method}")

    def percentage_of_clip_lengths(self, threshold: int) -> float:
        '''
        Checks what percentage of the CLIP data has a length above the specified threshold.

        Parameters
        ----------
        threshold: int
            The threshold, which corresponds to the nucleotide length to check for.

        Returns
        -------
        float
            Percentage of the data that is longer than the specified threshold.
        '''
        for method in self.unique_methods:
            lengths = [entry.clip_end - entry.clip_start for entry in self.filter_entries(method=method)]
            count_above = sum(1 for x in lengths if x > threshold)
            return round((count_above / len(lengths)) * 100, 3)

    def entries_of_software_for_rbp(self) -> dict[str, dict[str, int]]:
        """
        Creates a dictionary of all RBPs and the software tools used to determin the RBP binding sites.

        Returns
        -------
        dict[str, dict[str, int]]
            Dictionary with the number of entries for each unique software tool of every RBP.
        """
        rbps_tools_dict = {}
        for rbp in self.unique_rbp:
            clip_data = [clip for clip in self.clip_data if clip.rbp_name.upper() == rbp]
            tool_counts = Counter(clip.software for clip in clip_data)
            rbps_tools_dict[rbp] = dict(tool_counts)
        return rbps_tools_dict

    def filter_upper_decile_score(self) -> list[FullClipEntry]:
        """
        Filters the data with confidence scores in the top 10% for each experiment and software tool.

        Returns
        -------
        list[FullClipEntry]
            List of the filtered CLIP data.
        """
        filtered_clip_data = []
        for experiment in self.unique_accession_experiments:
            experiment_data = self._extract_upper_decile_scores_tool(experiment)
            filtered_clip_data.extend(experiment_data)
        return filtered_clip_data

    def _extract_upper_decile_scores_tool(self, experiment: str) -> list[FullClipEntry]:
        """
        Extracts the data with confidence scores in the top 10% of the specified experiment for all software tools.

        Parameters
        ----------
        experiment: str
            String of the experiment from the accession data.

        Returns
        -------
        list[FullClipEntry]
            List of the filtered data for the experiment.
        """
        filtered_data = []
        experiment_data: list[FullClipEntry] = [data for data in self.clip_data if data.accession_experiment == experiment]
        unique_software_tools = {data.software for data in experiment_data}
        for software_tool in unique_software_tools:
            filtered_data.extend(self._get_tool_specific_data(software_tool, experiment_data))
        return filtered_data
    
    def _get_tool_specific_data(self, software_tool: str, experiment_data: list[FullClipEntry]) -> list[FullClipEntry]:
        """
        Gets the data with confidence scores in the top 10% for the specified software tool.

        Parameters
        ----------
        software_tool: str
            Name of the software tool.
        experiment_data: list[FullClipEntry]
            List of the unfiltered experiment data.

        Returns
        -------
        list[FullClipEntry]
            Filtered data of the software tool.
        """
        tool_specific_data = [data for data in experiment_data if data.software == software_tool]
        tool_specific_data = sorted(tool_specific_data, key=lambda x: x.confidence_score, reverse=True)
        decile = int(len(tool_specific_data) * 0.1) -1
        if software_tool in self.software_tools:
            cutoff = tool_specific_data[decile].confidence_score
            decile_data = [data for data in tool_specific_data if data.confidence_score >= cutoff]
        else:
            cutoff = tool_specific_data[-decile].confidence_score
            decile_data = [data for data in tool_specific_data if data.confidence_score >= cutoff]
        return decile_data
                
    def _extract_clip_data(self, organism: str) -> list[FullClipEntry]:
        """
        Extracts the CLIP data from the .txt file of the crresponding organism.

        Parameters
        ----------
        organism: str
            Name string of the organism.

        Returns
        -------
        list[FullClipEntry]
            List of FullClipEntry data.
        """
        clip_data, gff3_file = [], f"data/datasets/{organism}/{organism}_clip.txt"
        with open(gff3_file, "r") as file:
            for line in file:
                entry = line.strip().split("\t")
                clip_entry = self._parse_clip_entry(entry)
                clip_data.append(clip_entry)
        return clip_data

    @staticmethod
    def _parse_clip_entry(entry: list[str]) -> FullClipEntry:
        """
        Parses an entry from a .txt file and turns it into a FullCLipEntry.

        Parameters
        ----------
        entry: list[str]
            Entry with the information as a list of strings.

        Returns
        -------
        FullClipEntry
            The entry as a FullClipEntry object.
        """
        method_software = entry[6].split(",")
        method, software = (method_software if len(method_software) == 2 
                            else (method_software[0], method_software[0]))
        accession = entry[8].split(",")
        accession_data, accession_experiment = (accession if len(accession) == 2 
                                                else (accession[0], ""))
        return FullClipEntry(
            clip_id="",
            start=int(entry[1]),
            end=int(entry[2]),
            strand_orientation=entry[4],
            rbp_name=entry[5],
            method=method,
            software=software,
            sample=entry[7],
            accession_data=accession_data,
            confidence_score=float(entry[9]),
            chromosome=entry[0],
            peak_id=entry[3],
            accession_experiment=accession_experiment,
            feature_types=[],
            sequence="",
            sequence_start=None,
            sequence_end=None
        )
    
class ClipSpeciesAnalyzer:
    """
    Class analyzing all species specified in the config file.

    Attributes
    ----------
    species: dict[str, str]
        Dictionary containing the species names.
    software_tools: dict[str, bool]
        Dictionary of software tools.
    """

    def __init__(self, config_file: str):
        """
        Initializes a ClipSpeciesAnalyzer object.

        Parameters
        ----------
        config_file: str
            Yaml file containing the species names.
        """
        with open(config_file, 'r') as file:
            config = yaml.safe_load(file)
            self.species: dict[str, str] = config["organisms"]

    def analyse_all_organisms(self) -> None:
        """
        Analyzes all organisms and prints the values, which are of interest.
        """
        # Use this method to analyze specific properties of the datasets.
        for species_name, species in self.species.items():
            print(f"Analyzing {species_name}:")
            data_analyzer = ClipDataAnalyzer(species, HIGH_CONFIDENCE_SCORER)
            #print(f"unique: chromosomes: {data_analyzer.unique_chromosomes}")
            #print(
            #    f"Unique chromosomes: {data_analyzer.unique_chromosomes}",
            #    f"Unique methods: {data_analyzer.unique_methods}",
            #    f"Unique software tools: {data_analyzer.unique_software_tools}",
            #    f"Unique RBPs: {data_analyzer.unique_rbp}",
            #    f"Number of unique accession experiments/samples: {len(data_analyzer.unique_accession_experiments)}",
            #    f"Number of entries: {len(data_analyzer.clip_data)}",
            #    f"Number of filtered: {len(data_analyzer.filter_upper_decile_score())}",
            #    f"Size of list object:", sys.getsizeof(data_analyzer), "bytes",
            #    sep="\n", end="\n\n"
            #    )
            #print(data_analyzer.entries_of_software_for_rbp(), "\n")
            #data_analyzer.all_data_clip_length_histogram()
            print(f"Percentage above 50 nt: {data_analyzer.percentage_of_clip_lengths(50)}%")