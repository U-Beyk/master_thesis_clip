'''
run_processing.py
=================

Module containg the function to run the processing stage of the pipeline.

author: U.B.
'''

from .fasta_builder import FastaBuilder

def run_processing() -> None:
    '''
    Runs the processing stage of the CLIP-seq pipeline 
    by building the FASTAs.
    '''
    builder = FastaBuilder()
    builder.build_fasta_organisms()