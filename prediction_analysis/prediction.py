import re
from dataclasses import dataclass, field

# TO-DO: Implement recognition of ambiguous motifs instead of hardcoding
AMBIGUOUS_MOTIFS = {"u": "GU", "g": "GT", "t": "GT"}

@dataclass(slots=False)
class CLIPmotiFoldPrediction:
    '''
    Standard blueprint for RNAmotiFold prediction.
    '''
    prediction_class: str
    free_energy: float
    mot_bracket: str
    rel_start_end: tuple[int, int]
    distance_to_mfe: float | None = field(init=False, default=None)

    @property
    def motifs_set(self) -> set[str]:
        '''Gets and returns a set of unique motifs in the prediction, replacing ambiguous lowercase motifs.'''
        motifs = set()
        for motif in self.prediction_class:
            if motif.islower() and motif in AMBIGUOUS_MOTIFS:
                motifs.update(AMBIGUOUS_MOTIFS[motif])
            else:
                motifs.update(motif)
        return motifs
    
    @property
    def mfe_distance_motif(self) -> dict[str, float]:
        if not self.motifs_set:
            return {"No motif": self.distance_to_mfe}
        return {motif: self.distance_to_mfe for motif in self.motifs_set}

    # TO-DO: Check indexing, if correct.
    @property
    def clip_overlap_motifs_set(self) -> None:
        clip_structure = self.mot_bracket[self.rel_start_end[0]: self.rel_start_end[1] + 1]
        motifs = set()
        for motif in clip_structure:
            if motif in ".()":
                continue
            if motif.islower() and motif in AMBIGUOUS_MOTIFS:
                motifs.update(AMBIGUOUS_MOTIFS[motif])
            else:
                motifs.update(motif)
        return motifs

    @property
    def mfe_distance_clip_motif(self) -> dict[str, float]:
        '''Assigns to each motif the value for the distance to the mfe.
        If there are no motifs, it returns a corresponding dictionary.'''
        if not self.clip_overlap_motifs_set:
            return {"No motif": self.distance_to_mfe}
        return {
            motif: self.distance_to_mfe 
            for motif in self.clip_overlap_motifs_set
        }
    
@dataclass(slots=False)
class CLIPmotiCesPrediction(CLIPmotiFoldPrediction):
    potential_motif_sequences: list[str] = field(init=False, default_factory=list)

    # TO-DO: Change or implement new method, that gets only the potential motifs, that overlap with CLIP data.
    def compute_potential_motifs(self, sequence) -> None:
        '''Computes the potential motif for the corresponding sequence of that prediction.'''
        for position in self.positions_without_motif:
            left_nt_pos = self._find_left_bracket_position(round(float(position)))
            right_nt_pos = self._find_right_bracket_position(round(float(position)))
            self.potential_motif_sequences.append(sequence[left_nt_pos: right_nt_pos])

    def _find_left_bracket_position(self, start_pos: int) -> int | None:
        '''Finds the nearest '(' bracket position to the left.'''
        for pos in range(start_pos, 0, -1):
            if self.mot_bracket[pos] == "(":
                return pos + 1
        return None

    def _find_right_bracket_position(self, start_pos: int) -> int | None:
        '''Finds the nearest ')' bracket position to the right.'''
        for pos in range(start_pos, len(self.mot_bracket)):
            if self.mot_bracket[pos] == ")":
                return pos
        return None

    @property
    def motifs_set(self) -> set[str]:
        '''Gets and returns a set of unique motifs, replacing ambiguous ones.'''
        motifs = set()
        for motif in str(self.prediction_class):
            if motif.isalpha():
                if motif.islower() and motif in AMBIGUOUS_MOTIFS:
                    motifs.update(AMBIGUOUS_MOTIFS[motif])
                else:
                    motifs.update(motif)
        return motifs
    
    @property
    def positions_without_motif(self) -> list[float]:
        '''finds and returns all the secondary structures without a motif.'''
        return re.findall(r'\d+\.?\d+(?=_)', self.prediction_class)