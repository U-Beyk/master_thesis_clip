'''
clip_processor.py
=================

Module with the class processing the CLIP data.

author: U.B.
'''

from collections.abc import Iterator

from .clip_entry import ClipEntry
from .clip_filter import ClipFilter
from ..utils import normalize_chr

class ClipProcessor:
    '''
    Handles reading and processing of CLIP data from a .txt file.

    Attributes
    ----------
    clip_file: str
        Path to the CLIP input file.
    _filter: ClipFilter
        Filtering logic for entry validation.
    '''

    def __init__(self, clip_file: str):
        '''
        Initiates a ClipProcessor object.

        Parameters
        ----------
        clip_file: str
            Path of the file with the CLIP data.
        '''
        self.clip_file = clip_file
        self._filter = ClipFilter(self._read_file())

    def iterate_clips(self) -> Iterator[ClipEntry]:
        '''
        Iterates over parsed clip entries that pass the configured filters.

        Yields
        ------
        ClipEntry
            A parsed and validated clip entry that passes all filter checks.
        '''
        for entry in self._read_file():
            parsed = self._parse_entry(entry)
            if not self._filter.passes(parsed):
                continue
            yield parsed

    def _read_file(self) -> Iterator[list[str]]:
        '''
        Reads the clip data file line by line and yields each entry as a list of fields.

        Yields
        ------
        list[str]
            A list of tab-separated values representing one raw clip record.
        '''
        with open(self.clip_file) as f:
            for line in f:
                yield line.strip().split("\t")

    def _parse_entry(self, entry: list[str]) -> ClipEntry:
        '''
        Parses a CLIP entry into a structured "ClipEntry" object.

        Parameters
        ----------
        entry: list[str]
            A single tab-split record representing one clip entry.

        Returns
        -------
        ClipEntry
            A "ClipEntry" object.
        '''
        chromosome = normalize_chr(entry[0])
        start, end = int(entry[1]), int(entry[2])
        strand, rbp_name = entry[4], entry[5]
        method_parts = entry[6].split(",")
        method, software = method_parts[0], method_parts[-1]
        sample = entry[7]
        acc_parts = entry[8].split(",")
        accession_data, accession_experiment = acc_parts[0], acc_parts[-1]
        score = float(entry[9])

        return ClipEntry(
                clip_id="",
                chromosome=chromosome,
                clip_start=start,
                clip_end=end,
                strand_orientation=strand,
                rbp_name=rbp_name,
                method=method,
                software=software,
                sample=sample,
                accession_data=accession_data,
                accession_experiment=accession_experiment,
                confidence_score=score,
                feature_types=[],
                sequence="",
                sequence_start=None,
                sequence_end=None,
                shuffled_sequence=""
            )