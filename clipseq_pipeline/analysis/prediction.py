import re
from dataclasses import dataclass, field

from ..constants import AMBIGUOUS_MOTIFS
from ..constants import CLIP_NT_WINDOW

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

    # TODO: Check indexing, if correct.
    # TODO: Refactor code. 
    @property
    def clip_overlap_motifs_set(self) -> None:
        clip_start, clip_end = self.rel_start_end
        clip_length = clip_end - clip_start + 1
        target_length = CLIP_NT_WINDOW

        # Extends evenly on both sides
        extend_total = max(0, target_length - clip_length)
        extend_left = extend_total // 2
        extend_right = extend_total - extend_left

        # Applies left extension, respecting lower bound (0)
        seq_start = max(0, clip_start - extend_left)
        # Applies right extension, respecting upper bound
        seq_end = min(len(self.mot_bracket), clip_end + extend_right)
        
        actual_right_ext = seq_end - clip_end  # actual extension on right side
        actual_left_ext = clip_start - seq_start  # actual extension on left side

        # Extends the other side when another side hit s a boundary
        if (actual_left_ext + actual_right_ext) < extend_total:
            remaining = extend_total - (actual_left_ext + actual_right_ext)

            # Extends the right side if possible
            if seq_end < len(self.mot_bracket):
                add_right = min(remaining, len(self.mot_bracket) - seq_end)
                seq_end += add_right
                remaining -= add_right

            # Extends the left side if possible
            if remaining > 0 and seq_start > 0:
                seq_start = max(0, seq_start - remaining)

        clip_structure = self.mot_bracket[seq_start:seq_end]
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
    
# TODO: Change or implement new method, that gets only the potential motifs, that overlap with CLIP data.
@dataclass(slots=False)
class CLIPmotiCesPrediction(CLIPmotiFoldPrediction):
    potential_motif_sequences: list[str] = field(init=False, default_factory=list)

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