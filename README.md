# Master thesis

**THIS IS A WORK IN PROGRESS!**

The datasets used in this master's thesis were downloaded on July 22, 2025.
We downloaded the CLIP data from the POSTAR3 database by requesting a bulk download:
http://postar.ncrnalab.org

The different genome and gene annotation files were downloaded from Ensembl.

For Homo sapiens (human; release: 114):
Genome: https://ftp.ensembl.org/pub/release-114/fasta/homo_sapiens/dna/
Gene annotations: https://ftp.ensembl.org/pub/release-114/gff3/homo_sapiens/

For Mus musculus (mouse; release: 114):
Genome: https://ftp.ensembl.org/pub/release-114/fasta/mus_musculus/dna/
Gene annotations: https://ftp.ensembl.org/pub/release-114/gff3/mus_musculus/

For Danio rerio (zebrafish; release: 114):
Genome: https://ftp.ensembl.org/pub/release-114/fasta/danio_rerio/dna/
Gene annotations: https://ftp.ensembl.org/pub/release-114/gff3/danio_rerio/

For Drosophila Melanogaster (fly; release: 114):
Genome: https://ftp.ensembl.org/pub/release-114/fasta/drosophila_melanogaster/dna/
Gene annotations: https://ftp.ensembl.org/pub/release-114/gff3/drosophila_melanogaster/

For Caenorhabditis elegans (worm; release: 114):
Genome: https://ftp.ensembl.org/pub/release-114/fasta/caenorhabditis_elegans/dna/
Gene annotations: https://ftp.ensembl.org/pub/release-114/gff3/caenorhabditis_elegans/

For Saccharomyces cerevisiae (yeast; release: 114):
Genome: https://ftp.ensembl.org/pub/release-114/fasta/saccharomyces_cerevisiae/dna/
Gene annotations: https://ftp.ensembl.org/pub/release-114/gff3/saccharomyces_cerevisiae/

For Arabidopsis thaliana (release: 61):
Genome: https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-61/fasta/arabidopsis_thaliana/dna/
Gene annotations: https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-61/gff3/arabidopsis_thaliana/

The genome files containing all the DNA data for each organism were downloaded and unpacked, with file types left unchanged after extraction. 
A dedicated folder was created for each organism within the directory data/datasets/.
For example, the folder for Homo sapiens was named files/dataset/homo_sapiens, and all relevant files for that organism were placed inside this folder.
To maintain consistency, all files were renamed using a standardized naming convention based on the folder name. 
The following format was used: {folder_name}_genome, {folder_name}_annotations, {folder_name}_clip (for the CLIP data)
For instance, in the Homo sapiens folder, the files were named as follows:
homo_sapiens_genome, homo_sapiens_annotations and homo_sapiens_clip.

A config.yaml file was created to include these organisms and other importatn settings.

Before running the project, the required packages should be installed via the requirements.txt with following command:
**pip install -r requirements.txt**

The analysis can be run with: **python3 main.py**

The tool can use up a lot of RAM, so having at around ~30GB free RAM is recommended,
otherwise the number of workers in the config.yaml can be decreased. 
If the amount of free RAM is much higher the number can also be increased.

To be on the safe site, use following command when running the analysis:
**systemd-run --scope -p MemoryMax=28G --user python3 main.py**
With this command, it is possible to specify a max memory usage. If that number is exceeded,
the program will be terminated.