'''
Contains all classes with the dataframes to be analyzed and its methods to filter necessary information.

author: U.B.
'''

from dataclasses import dataclass

from prediction_analysis.dataframe_assembler import DataframeAssembler
from prediction_analysis.records_assembler import CLIPmotiFoldRecordsAssembler, CLIPmotiCesRecordsAssembler

@dataclass
class AnalyzerConfig:
    '''Dataclass containing the information of the fasta file and the RNAmotiFold and RNAHeliCes/RNAmotiCes prediction files.'''
    fasta_path: str
    motifold_csv_path: str

class Analyzer:
    '''Generic class for the different RNA assemblies.'''

    def __init__(self, config: AnalyzerConfig) -> None:
        '''Initializes an Analyzer object with an RNAmotiFold assembly and an RNAHeliCes/RNAmotiCes assembly.'''
        self.rna_motifold = CLIPmotiFoldRecordsAssembler(
            DataframeAssembler(config.fasta_path, config.motifold_csv_path)
        )

    def filter_motifold_mfe_range(self) -> CLIPmotiFoldRecordsAssembler:
        '''Filters the predictions with a mfe range in the RNAmotiFold assembly.'''
        return self.rna_motifold.filter_records(mfe_range=True)
    
    def filter_motifold_by_rbp(self, rbp: str, mfe_range: bool = False) -> CLIPmotiFoldRecordsAssembler:
        '''Filters the RNAmotiFold records corresponding to a specified CRISPR subtype.
        Also allows to filter predictions within a mfe range.'''
        return self.rna_motifold.filter_records(
            condition=lambda record: rbp == record.rbp_name,
            mfe_range=mfe_range
        )

    def __repr__(self) -> str:
        '''Represents the object and its details as a string.'''
        return (f"{self.__class__.__name__}("
                f"rna_motifold_assembly={self.rna_motifold}")