# All organisms to include in the analysis of the CLIP-seq data, some of them csan be commented out, to test stuff
ORGANISMS_TO_EXAMINE = {
    #"Arabidopsis thaliana": "arabidopsis_thaliana",
    #"Caenorhabditis elegans": "caenorhabditis_elegans",
    #"Danio rerio": "danio_rerio",
    "Drosophila melanogaster": "drosophila_melanogaster",
    #"Homo sapiens": "homo_sapiens",
    #"Mus musculus": "mus_musculus",
    #"Saccharomyces cerevisiae": "saccharomyces_cerevisiae"
}

# Length of the CLIP data to include. Everything above will be excluded from the processing.
CLIP_NT_LENGTH = 50

# Determines the length of the final RNA sequence to predict
RBP_NT_LENGTH = 150

# These are peak-calling software tools that score CLIP-seq data where a higher score equal a higher confidence.
HIGH_CONFIDENCE_SCORER = [
    "eCLIP", "PARalyzer", "PureCLIP", "CTK", "CIMS", "MiClip",
    "Piranha_0.01", "PIP-seq"
]


# Maps ambiguous chromosome names to normalize them.
AMBIGUOUS_CHR_MAPPING = {
    # Mitochondrial variants
    "M": "MT",
    "MTDNA": "MT",
    "MITO": "MT"
}

# Correpsonds to the number of workers for the RNA predictions.
RNA_PREDICTION_WORKER = 12

# This value determines the number of the lowest free energy values of the classes
RNA_PREDICTION_K_VALUE = 10


# Ambiguous motifs from the RNAmotiFold prediction, where lower case means that the motif can equal to two different motifs.
AMBIGUOUS_MOTIFS = {"u": "GU", "g": "GT", "t": "GT"}

# The max CLIP window in the RNA sequence to use for analysing
CLIP_NT_WINDOW = 50