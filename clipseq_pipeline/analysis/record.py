import copy
from typing import Self
from dataclasses import dataclass, field

from .prediction import CLIPmotiFoldPrediction, CLIPmotiCesPrediction 

@dataclass(slots=False)
class CLIPmotiFoldRecord:
    sequence_id: str
    clip_start: int
    clip_end: int
    rbp_name: str
    sequence_start: int
    sequence_end: int
    feature_types: list[str]
    sequence: str
    predictions: list[CLIPmotiFoldPrediction] = field(default_factory=list)

    @classmethod
    def from_header(cls, sequence_header: str, sequence: str) -> Self:
        """Parses a sequence header and returns a CLIPRecord instance."""
        id_part = sequence_header.split("|")
        sequence_id = id_part[0].split("id:")[1]
        clip_start, clip_end = map(int, id_part[1].split("clip_range:")[1].split("-"))
        rbp_name = id_part[2].split("rbp_name:")[1]
        sequence_start, sequence_end = map(int, id_part[3].split("seq_range:")[1].split("-"))
        feature_types = id_part[4].split("features:")[1].split(",")
        return cls(sequence_id, clip_start, clip_end, rbp_name, sequence_start, sequence_end, feature_types, sequence)
    
    @property
    def rel_clip_start_end(self) -> tuple[int, int]:
        """Returns the clip start and end positions relative to the sequence start."""
        return (self.clip_start - self.sequence_start,
                self.clip_end - self.sequence_start)

    @property
    def mfe_prediction(self) -> CLIPmotiFoldPrediction | None:
        '''Calculates and returns the prediction with the mfe value.'''
        return min(
            self.predictions,
            key=lambda prediction: prediction.free_energy,
            default=None
            )
    
    @property
    def mfe_value(self) -> float | None:
        '''Returns the mfe value as a float.'''
        return self.mfe_prediction.free_energy if self.mfe_prediction else None
    
    @property
    def prediction_number(self) -> int:
        '''Returns the total number of all predictions.'''
        return len(self.predictions)
    
    @property
    def motifs_set(self) -> set[str]:
        motifs = {motif for p in self.predictions for motif in p.motifs_set}
        return motifs or {"No motif"}
    
    @property
    def mfe_motif_distances(self) -> list[dict[str, float]]:
        '''Returns a list of all motifs '''
        return [
            prediction.mfe_distance_motif
            for prediction in self.predictions
            ]
    
    @property
    def clip_overlap_motifs_set(self) -> set[str]:
        motifs = {motif for p in self.predictions for motif in p.clip_overlap_motifs_set}
        return motifs or {"No motif"}
    
    @property
    def mfe_clip_motif_distances(self) -> list[dict[str, float]]:
        return [
            prediction.mfe_distance_clip_motif
            for prediction in self.predictions
        ]
    
    def filter_predictions(self, mfe_range: float) -> Self:
        '''Filters predictions within a specified mfe range from the mfe and returns a new RNARecord instance.'''
        mfe = self.mfe_value
        new_record = copy.deepcopy(self)
        new_record.predictions = [
            pred for pred in self.predictions if pred.free_energy <= mfe + abs(mfe_range)
        ]
        new_record._update_mfe_distance()
        return new_record
    
    def add_prediction(self, prediction: CLIPmotiFoldPrediction) -> None:
        '''Adds a prediction and updates its distances to the lowest mfe value.'''
        self.predictions.append(prediction)
        self._update_mfe_distance()

    def _update_mfe_distance(self) -> None:
        '''Updates and calculates the distance to the mfe value for all predictions.'''
        if not self.predictions:
            return
        min_mfe = self.mfe_value
        for prediction in self.predictions:
            prediction.distance_to_mfe = prediction.free_energy - min_mfe

# TODO: Change or implement new method, that gets only the potential motifs, that overlap with CLIP data.
@dataclass(slots=False)
class CLIPmotiCesRecord(CLIPmotiFoldRecord):
    predictions: list[CLIPmotiCesPrediction] = field(default_factory=list) 

    @classmethod
    def from_header(cls, sequence_header: str, sequence: str) -> Self:
        """Parses a sequence header and returns a CLIPmotiCesRecord instance."""
        id_part = sequence_header.split("|")
        sequence_id = id_part[0].split("id:")[1]
        clip_start, clip_end = map(int, id_part[1].split("clip_range:")[1].split("-"))
        rbp_name = id_part[2].split("rbp_name:")[1]
        sequence_start, sequence_end = map(int, id_part[3].split("seq_range:")[1].split("-"))
        feature_types = id_part[4].split("features:")[1].split(",")
        return cls(sequence_id, clip_start, clip_end, rbp_name, sequence_start, sequence_end, feature_types, sequence)

    @property
    def mfe_prediction(self) -> CLIPmotiCesPrediction | None:
        '''Calculates and returns the prediction with the mfe value.'''
        return min(
            self.predictions,
            key=lambda prediction: prediction.free_energy,
            default=None
            )

    @property
    def potential_motifs_set(self) -> set[str]:
        '''Returns all the potential motifs that were computed via positional abstraction.'''
        potential_motifs = set()
        for prediction in self.predictions:
            prediction.compute_potential_motifs(self.sequence)
            potential_motifs.update(prediction.potential_motif_sequences)
        return potential_motifs

    def add_prediction(self, prediction: CLIPmotiCesPrediction) -> None:
        '''Adds a prediction and updates its distances to the lowest mfe value.'''
        self.predictions.append(prediction)
        self._update_mfe_distance()