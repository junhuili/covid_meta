# covid_meta

# System Requirements
TrimGalore v0.6.7, HUMAnN 3, MetaPhlAn 3 and Kraken2 can be installed in Conda; SPINGO v1.3 in Linux system.

# Softwares
QC with TrimGalore v0.6.7 | Software requirements: Cutadapt and FastQC | Installation guide: https://github.com/FelixKrueger/TrimGalore | Instructions for use: https://github.com/FelixKrueger/TrimGalore/blob/master/Docs/Trim_Galore_User_Guide.md

16S taxonomy using SPINGO v1.3 | Software requirements, Installation guide, Demo, and Instructions for use: https://github.com/GuyAllard/SPINGO 

Shotgun function using HUMAnN 3 | Software requirements, Installation guide, Demo, and Instructions for use: https://github.com/biobakery/humann 

Shotgun taxonomy using MetaPhlAn 3 | Software requirements, Installation guide, Demo, and Instructions for use: https://github.com/biobakery/MetaPhlAn 

Shotgun taxonomy using Kraken2 | Software requirements, Installation guide, Demo, and Instructions for use: https://github.com/DerrickWood/kraken2/wiki/Manual

R v3.6.2 | Installation guide, Demo, and Instructions for use: https://cran.r-project.org/manuals.html 

# Databases
Human and PhiX genome database: https://doi.org/10.6084/m9.figshare.12525242.v1

HUMAnN3: http://cmprod1.cibio.unitn.it/databases/HUMAnN/full_mapping_v201901.tar.gz, http://cmprod1.cibio.unitn.it/databases/HUMAnN/uniref90_annotated_v201901b.tar.gz

MetaPhlAn3: 10.5281/zenodo.3955744

Kraken2: bacteria, archaea, fungi, viral, protozoa dowloaded from NCBI on 2022/02/18 (see kraken2.sh in bioinformatics_statistics)

SPINGO: database/RDP_11.2.species.zip, taxonomy.zip


# Microbiome test data
The example microbiome data PRJNA740067 (shotgun metagenome) and PRJNA684070 (16S rRNA gene amplicon) can be downloaded from NCBI or EBI.

# Workflow
1. Bioinformatic anlysis pipeline (details for softwares and parameters) in the bioinformatics_statistics folder: bioinformatics.sh, kraken2.sh.

Shotgun metagenomes were quality-filtered using TrimGalore and further filtered to remove human and PhiX sequences using Kraken2. The decontaminated sequences were reprocessed using MetaPhlAn3 for species-level taxonomic profiling and HUMAnN3 for functional profiling that were further regrouped to KEGG Orthologs (KOs), Enzyme Commission number (EC number), MetaCyc pathway (Pathway), and carbohydrate active enzymes (CAZymes). Besides, the decontaminated sequences were reprocessed using Kraken2 for species-level taxonomic profiling of 5 microbial kingdoms against the constructed database using all NCBI genomes from five kingdoms (i.e., bacteria, archaea, fungi, viruses, and parasites).

16S sequences were also quality-filtered using TrimGalore and then reprocessed using SPINGO32 for taxonomic profiling.

2. Statistics and visualization: 

R codes for data analysis and visualisation in the bioinformatics_statistics folder: NMDS.R, PERMANOVA.R, dysbiosis.R, association_heatmap.R, RandomForest_COVIDSeverityIndex.R

R Package Dependencies: ape, ggplot2, vegan, BiodiversityR, MASS, pROC, randomForest, ggplot2

