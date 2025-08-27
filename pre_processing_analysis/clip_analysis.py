"""
clip_analysis.py
~~~~~~~~~~~~~~~~

This module contains all classes for analyzing the CLIP data before processing.

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

from processing.clip_annotation import ClipEntry

@dataclass(frozen=True,slots=True)
class FullClipEntry(ClipEntry):
    """
    FullClipEntry class representing the blueprint for clip data entries.
    This class includes atributes like chromosome, peak_id and accession_experiment.
    The rest is inherited from the ClipEntry class.

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
    chromosome: str
        Number/name of the chromosome.
    peak_id: str
        ID of the peak in te POSTAR3 data.
    acession_experiment: str
        Experiment accession data string.
    """
    chromosome: str
    peak_id: str
    accession_experiment: str

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
            "start": self.start,
            "end": self.end,
            "peak ID": self.peak_id,
            "strand orientation": self.strand_orientation,
            "RNA binding protein": self.rbp_name,
            "method": self.method,
            "software": self.software,
            "sample/tissue": self.sample,
            "accession of raw data": self.accession_data,
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
        Dictionary of software tools.
        True means that a higher score
        corresponds to a higher confidence.
        False means that a lower score
        corresponds to a higher confidence.
    clip_data: list[FullClipEntry]
        List of CLIP data.
    """

    __slots__ = ("organism", "clip_data", "software_tools")
    def __init__(self, organism: str, software_tools: dict[str, bool]):
        """
        Initializes a ClipDataAnalyzer instance.

        Parameters
        ----------
        organism: str
            Organism the CLIP data to extract from.
        software_tools: dict[str, bool]
            Dictionary of software tools to include in the analysis.
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
        plt.savefig(f"plots/pre_processing_plots/{self.organism}/{file_name}.png")
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
    
    def all_data_histogram(self) -> None:
        """
        Creates histograms for each unique software tool.
        """
        for tool in self.unique_software_tools:
            score_data = [entry.confidence_score for entry in self.filter_entries(software=tool)]
            if score_data:
                self.save_histogram(score_data, tool)

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
        if self.software_tools.get(software_tool):
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
            accession_experiment=accession_experiment
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
            self.software_tools: dict[str, bool] = config["included_software_tools"]

    def analyse_all_organisms(self) -> None:
        """
        Analyzes all organisms and prints the values, which are of interest.
        """
        # Use this method to analyze specific properties of the datasets.
        for species_name, species in self.species.items():
            print(f"Analyzing {species_name}:")
            data_analyzer = ClipDataAnalyzer(species, self.software_tools)
            print(
                f"Unique chromosomes: {data_analyzer.unique_chromosomes}",
                f"Unique methods: {data_analyzer.unique_methods}",
                f"Unique software tools: {data_analyzer.unique_software_tools}",
                #f"Unique RBPs: {data_analyzer.unique_rbp}",
                f"Number of unique accession experiments/samples: {len(data_analyzer.unique_accession_experiments)}",
                f"Number of entries: {len(data_analyzer.clip_data)}",
                f"Number of filtered: {len(data_analyzer.filter_upper_decile_score())}",
                sep="\n", end="\n\n"
                )
            #print(data_analyzer.entries_of_software_for_rbp(), "\n")
            #data_analyzer.all_data_histogram()