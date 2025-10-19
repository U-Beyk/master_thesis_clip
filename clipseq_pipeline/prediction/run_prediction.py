'''
run_prediction.py
=================

Module containg the function to run the prediction stage of the pipeline.

author: U.B.
'''

from .rna_predicter import RNAPredicter

def run_prediction():
    '''
    Runs the prediction stage of the CLIP-seq pipeline 
    by building the FASTAs.
    '''
    predicter = RNAPredicter()
    predicter.run_predictions()