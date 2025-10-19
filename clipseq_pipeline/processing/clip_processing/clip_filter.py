'''
clip_filter.py
==============

This modules contains the ClipFilter class, which is responsible
for filtering valid CLIP entries to further process.

author: U.B.
'''

from collections import defaultdict
from collections.abc import Iterator

from ...constants import CLIP_NT_LENGTH, HIGH_CONFIDENCE_SCORER
from .clip_entry import ClipEntry


class ClipFilter:
    '''
    Class handling filtering logic.
    This includes the CLIP data having the right length,
    having a score value and the score value being in the
    cutoff range of that experiment and used software tool.

    Attributes
    ----------
    _high_confidence_scorer: list[str]
        List of the peak calling software, where a higher score
        equals to a higher confidence of the peak being a 
        RNA-binding-protein binding site.
    _cutoffs: dict[str, dict[str, float]]
        Dictionary with the cutoffs for all experiments and of each experiment
        for all peak calling software.
    '''

    def __init__(self, entries: Iterator[list[str]]):
        '''
        Instantiates a ClipFilter object.

        Parameters
        ----------
        entries: Iterator[list[str]]
            An iterator that yields lists of strings, where each list represents
            an entry from the CLIP dataset.
        '''
        self._high_confidence_scorer = HIGH_CONFIDENCE_SCORER
        self._cutoffs = self._compute_cutoffs(entries)

    def _compute_cutoffs(self, entries: Iterator[list[str]]) -> dict[str, dict[str, float]]:
        '''
        Compute cutoff values for each experiment/software combination:
        {"experiment1": {"software1": float, "software2": float}, "experiment2": {"software1": float}}

        Parameters
        ----------
        entries: Iterator[list[str]]
            An iterator that yields lists of strings, where each list represents
            an entry from the CLIP dataset.
        
        Returns
        -------
        dict[str, dict[str, float]]
            Dictionary with cutoff values for the experiment/software combinations.        
        '''
        scores = defaultdict(lambda: defaultdict(list))
        for entry in entries:
            experiment, software, score = self._parse_for_cutoff(entry)
            scores[experiment][software].append(score)

        cutoffs = {}
        for experiment, software_tool in scores.items():
            cutoffs[experiment] = {
                software: self._calculate_cutoff(values, software) 
                for software, values in software_tool.items()
            }
        return cutoffs

    @staticmethod
    def _parse_for_cutoff(entry: list[str]) -> tuple[str, str, float]:
        '''
        Parses the experiment, software and confidence score of a CLIP entry
        for further calculation of the cutoff.

        Parameters
        ----------
        entry: list[str]
            A list of the values of a CLIP entry.

        Returns
        -------
        tuple[str, str, float]
            A tuple containing the experiment string,
            the peak calling software string and
            the confidence score as a float.
        '''
        score = float(entry[9])
        acc_parts = entry[8].split(",")
        experiment = acc_parts[1] if len(acc_parts) == 2 else acc_parts[0]
        meth_parts = entry[6].split(",")
        software = meth_parts[1] if len(meth_parts) == 2 else meth_parts[0]
        return experiment, software, score

    def _calculate_cutoff(self, scores: list[float], software: str) -> float:
        '''
        Compute 10% cutoff for given score list and software.

        Parameters
        ----------
        scores: list[float]
            List of score values.
        software: str
            Peak calling software string.

        Returns
        -------
        float
            Calculated value for the cutoff of the corresponding software.
        '''
        scores.sort(reverse=True)
        index = max(1, int(len(scores) * 0.1)) - 1
        if software in self._high_confidence_scorer:
            # Use top 10% cutoff
            return scores[index]
        else:
            # Use bottom 10% cutoff
            # NOTE: scores[-idx] is wrong when idx=0, as -0 == 0 (top score), not bottom
            return scores[-(index + 1)]

    def passes(self, entry: ClipEntry) -> bool:
        '''
        Checks whether a given "ClipEntry" passes all configured filters.

        Parameters
        ----------
        entry: ClipEntry
            A parsed clip entry.

        Returns
        -------
        bool
            True if the entry meets the filtering criteria; False otherwise. 
        '''
        cutoff = self._cutoffs.get(entry.accession_experiment, {}).get(entry.software)
        if cutoff is None:
            return False
        
        if (entry.clip_end - entry.clip_start) > CLIP_NT_LENGTH:
            return False

        if entry.software in self._high_confidence_scorer:
            return entry.confidence_score >= cutoff
        else:
            return entry.confidence_score <= cutoff