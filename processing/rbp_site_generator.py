# TO-DO: Has to be adjusted, since grouping the CLIPs into RBP sites is not an option.

from dataclasses import dataclass
from collections import defaultdict

from processing.clip_processing import ClipEntry, ClipProcessing
from processing.genome_processing import GenomeProcessing
from processing.gff3_processing import Gff3Processing

@dataclass(slots=True)
class RbpSite:
    chromosome: str
    start: int
    end: int
    strand_orientation: str
    clip_data: list[ClipEntry]
    feature_types: list[str]
    sequence: str

    def to_fasta(self) -> str:
        unique_rbps = sorted({clip.rbp_name for clip in self.clip_data})
        rbp_info = ",".join(unique_rbps)
        header = (
            f">chromosome:{self.chromosome}|start:{self.start}|end:{self.end}"
            f"|strand_orientation:{self.strand_orientation}|RBPs={rbp_info}"
        )
        return f"{header}\n{self.sequence}\n"
    
class RbpSiteGenerator:
    def __init__(self, organism: str):
        self.organism = organism
        self.clip_data = ClipProcessing(f"data/datasets/{self.organism}/{self.organism}_clip.txt")
        self.genome = GenomeProcessing(f"data/datasets/{self.organism}/{self.organism}_genome.fa")
        self.gff3_index = Gff3Processing(f"data/datasets/{self.organism}/{self.organism}_annotations.gff3")

    def _group_clips_by_chromosome_and_strand(self) -> dict[tuple[str, str], list[ClipEntry]]:
        grouped_clips = defaultdict(list)
        for clip in self.clip_data.iterate_clips():
            key = (clip.chromosome, clip.strand_orientation)
            grouped_clips[key].append(clip)
        return grouped_clips
    
    def _group_clips_into_sites(self, clips, chromosome: str, strand_orientation: str) -> list[RbpSite]:
        rbp_sites: list[RbpSite] = []
        for clip in clips:
            if not rbp_sites:
                rbp_sites.append(self._create_rbp_site(chromosome, strand_orientation, clip))
                continue
            last_rbp_site = rbp_sites[-1]
            if self._overlaps(last_rbp_site, clip):
                self._merge_clip_into_site (last_rbp_site, clip)
            else:
                rbp_sites.append(self._create_rbp_site(chromosome, strand_orientation, clip))
        return rbp_sites

    @staticmethod
    def _create_rbp_site(chromosome: str, strand_orientation: str, clip: ClipEntry) -> RbpSite:
        return RbpSite(
            chromosome=chromosome,
            start=clip.start,
            end=clip.end,
            strand_orientation=strand_orientation,
            clip_data=[clip],
            feature_types=[],
            sequence=""
        )

    @staticmethod
    def _overlaps(site: RbpSite, clip: ClipEntry) -> bool:
        return (
            (clip.start >= site.start and clip.end <= site.end) or
            (site.start >= clip.start and site.end <= clip.end)
        )

    @staticmethod
    def _merge_clip_into_site (site: RbpSite, clip: ClipEntry) -> None:
        site.start = min(site.start, clip.start)
        site.end = max(site.end, clip.end)
        site.clip_data.append(clip)

    def _assign_features(self, site: RbpSite) -> bool:
        features = []
        if self.gff3_index:
            features = self.gff3_index.get_features(
                site.chromosome, site.strand_orientation, site.start, site.end
            )
        if not features:
            return False
        site.feature_types = features
        return True

    def _assign_sequence(self, site: RbpSite) -> bool:
        if not self.genome:
            return False
        seq_start = max(1, site.start - 50)
        seq_end = site.end + 50
        site.sequence = self.genome.get_sequence(
            site.chromosome, seq_start, seq_end, site.strand_orientation
        ) or ""
        return bool(site.sequence)

    def iterate_rbpsites(self):
        grouped_clips = self._group_clips_by_chromosome_and_strand()
        for (chromosome, strand_orientation), clips in grouped_clips.items():
            clips.sort(key=lambda clip: (clip.start, -clip.end))
            rbp_sites = self._group_clips_into_sites(clips, chromosome, strand_orientation)
            for site in rbp_sites:
                if not self._assign_features(site):
                    continue
                if not self._assign_sequence(site):
                    continue
                yield site

    def write_fasta(self) -> None:
        with open(f"data/fasta_files/{self.organism}_rbp_sites.fasta", "w") as file:
            for rbp_site in self.iterate_rbpsites():
                rbp_fasta_string = rbp_site.to_fasta()
                file.write(rbp_fasta_string)