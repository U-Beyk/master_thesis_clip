"""
data_mapper.py
~~~~~~~~~~~~~~

Contains all classes to map the CLIP data and sequences to specified transcripts.

Classes
-------
TranscriptProcessor
    Maps the CLIP data and sequence to their correpsonding transcript of a chromosome.
    Outputs the mapped data of a chromosome in a .json file.
PipelineRunner
    Pipeline for the mapping of the data for the specified organism.
SpeciesProcessor
    Maps all data for all chromosome of all specified organisms in the config.yaml file.

author: U.B.
"""

import os
import yaml
import orjson
from concurrent.futures import ProcessPoolExecutor, as_completed

from processing.genome_processing import Genome
from processing.clip_annotation import ClipParser
from processing.transcript_annotation import TranscriptAnnotationParser, TranscriptAnnotation

class TranscriptProcessor:
    """
    Processes and maps the sequence and CLIP data to their corresponding transcripts of the specified chromosome.

    Attributes
    ----------
    gff3_file: str
        The path to the GFF3 file with the transcript annotations.
    clip_file: str
        The path to the .txt CLIP file with the CLIP data.
    genome_fasta: str
        The path to the fasta file with the genome of the organism.
    output_folder: str
        The path of the folder the .json files are saved in.
    included_features: str
        The feature types to be included. Features not in this list are excluded and not parsed.
    """
    def __init__(self, organism: str, included_features: list[str]):
        """
        Initializes a TranscriptProcessor object to parse and map the CLIP data
        and the sequence to its corresponding transcripts of a specified organism.

        Parameters
        ----------
        organism: str
            The specified organism string.
        included_features: str
            All features to be included in the parsing and mapping of transcripts.
        """
        self.gff3_file = f"data/datasets/{organism}/{organism}_annotations.gff3"
        self.clip_file = f"data/datasets/{organism}/{organism}_clip.txt"
        self.genome_fasta = f"data/datasets/{organism}/{organism}_genome.fa"
        self.output_folder = f"data/datasets/{organism}/processed_chromosomes"
        self.included_features = included_features

    def _map_clip_to_transcripts(self, chromosome: str) -> dict[str, TranscriptAnnotation]:
        """
        Maps the CLIP data to its corresponding transcripts they overlap with of a chromosome.

        Parameters
        ----------
        chromosome: str
            The chromosome the transcripts and CLIP data are located in.

        Returns
        -------
        dict[str, TranscriptAnnotation]
            Dictionary with the transcript IDs as its keys and the annotated transcripts
            with their CLIP data as the values.
        """
        parser = TranscriptAnnotationParser(self.gff3_file, self.included_features)
        transcripts, tree = parser.build_interval_tree_for_chrom(chromosome)
        clip_data = ClipParser(self.clip_file)
        clips = clip_data.filter_clips_of_chromosome(chromosome)
        for clip in clips:
            overlapping = tree.overlap(clip.start, clip.end)
            for iv in overlapping:
                transcript_id = iv.data
                transcripts[transcript_id].append_clip(clip)
        return transcripts
    
    def _set_sequence_of_transcripts(self, transcripts: dict[str, TranscriptAnnotation]) -> None:
        """
        Stores the sequence of the transcripts in their sequence attribute.

        Parameters
        ----------
        transcripts: dict[str, TranscriptAnnotation]
            Dictionary with annotated transcripts. Keys correspond to their IDs
            and the values correspond to the annotated transcript.
        """
        # Instantiate Genome inside worker to avoid pickling BufferedReader errors!
        genome = Genome(self.genome_fasta)
        for transcript in transcripts.values():
            transcript_sequence = genome.extract_sequence(
                transcript.chromosome,
                transcript.start,
                transcript.end,
                transcript.strand_orientation
            )
            transcript.sequence = transcript_sequence

    def process_chromosome(self, chromosome: str) -> tuple[str, int]:
        """
        Maps and processes all CLIP data and sequences to their corresponding transcripts of a specified chromosome.
        The mapped data is stored as a .json file.

        Parameters
        ----------
        chromosome: str
            The specified chromosome.

        Returns
        -------
        tuple[str, int]
            Tuple containing:
            str: The chromosome name/number.
            int: Number of transcripts.
        """
        transcripts = self._map_clip_to_transcripts(chromosome)
        self._set_sequence_of_transcripts(transcripts)
        filtered = {
            transcript_id: transcript.to_dict()
            for transcript_id, transcript in transcripts.items() 
            if transcript.clip_data and transcript.sequence
            }
        os.makedirs(self.output_folder, exist_ok=True)
        out_path = os.path.join(self.output_folder, f"{chromosome}_transcripts_with_clip.json")
        with open(out_path, "wb") as out_file:
            out_file.write(orjson.dumps(filtered))
        return chromosome, len(filtered)


class PipelineRunner:
    """
    Class of the Pipeline to run the mapping and processing of the transcripts for the specified organism.

    Attributes
    ----------
    organism: str
        Organism the data is from to process.
    included_features: list[str]
        All features to be included in the parsing and mapping.
    max_workers: int
        Number of processor cores and processes run at the same time for the mapping of the chromosomes.
    """
    def __init__(self, organism: str, included_features: list[str], max_workers: int):
        """
        Initializes a PipelineRunner instance.

        Parameters
        ----------
        organism: str
            Organism to process and map its transcripts.
        included_features: list[str]
            All features to be included in the parsing and mapping.
        max_workers: int
            Number of processor cores and processes run at the same time.
        """
        self.organism = organism
        self.included_features = included_features
        self.max_workers = max_workers

    def _get_valid_chromosomes(self) -> list[str]:
        """
        Filters out all valid chromosomes in the GFF3 file of the organism.

        Returns
        -------
        list[str]
            A list with all valid chromosomes.
        """
        parser = TranscriptAnnotationParser(
            f"data/datasets/{self.organism}/{self.organism}_annotations.gff3", 
            self.included_features
            )
        valid_chromosomes = parser.extract_chromosomes()
        print(f"Processing chromosomes: {sorted(valid_chromosomes)}")
        return valid_chromosomes

    def _process_chromosome_worker(self, chromosome: str) -> tuple[str, int]:
        """
        Initializes a TranscriptProcessor object and maps the data of a chromosome.
        These can be done in parallel by several workers.

        Parameters
        ----------
        chromosome: str
            The chromosome the data is mapped of.
        
        Returns
        -------
        tuple[str, int]
            Tuple containing:
            str: The chromosome.
            int: Number of transcripts.
        """
        processor = TranscriptProcessor(
            organism=self.organism,
            included_features=self.included_features,
        )
        return processor.process_chromosome(chromosome)

    def run(self) -> None:
        """
        Runs the pipeline by filtering valid chromosomes and mapping the data of these chromosomes.
        Several chromosomes can be processed at the same time.
        """
        valid_chromosomes = self._get_valid_chromosomes()
        with ProcessPoolExecutor(max_workers=self.max_workers) as executor:
            futures = {
                executor.submit(self._process_chromosome_worker, chromosome): chromosome
                for chromosome in valid_chromosomes
            }
            for future in as_completed(futures):
                chromosome = futures[future]
                try:
                    chrom, count = future.result()
                    print(f"Finished chromosome {chrom} with {count} transcripts with CLIP data.")
                except Exception as e:
                    print(f"Error processing chromosome {chromosome}: {e}")

class SpeciesProcessor:
    """
    Processes and maps the data of all species specified in the config file.

    Attributes
    ----------
    species: dict[str, str]
        Dictionary with all species to process.
    included_features: list[str]
        All features to be included in the parsing and mapping of all organisms.
    max_workers: int
        Number of processor cores and processes run at the same time.
    """
    def __init__(self, config_file: str):
        """
        Instantiates a SpeciesProcessor object processing all organisms.

        Parameters
        ----------
        config_file: str
            Path of the config file with all the settings.
        """
        with open(config_file, 'r') as file:
            config = yaml.safe_load(file)
        self.species: dict = config["organisms"]
        self.included_features = config["included_features"]
        self.max_workers = config["number_of_workers_per_organism"]

    def process_species(self) -> None:
        """
        Processes and runs the pipeline for every species.
        """
        for species_name, species in self.species.items():
            print(f"\nRunning pipeline for {species_name}")
            pipeline = PipelineRunner(species, self.included_features, self.max_workers)
            pipeline.run()