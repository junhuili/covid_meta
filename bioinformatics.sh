##########################################
# 16S amplicon
##########################################

##############
# single read
##############

# QC
####
trim_galore --cores 8 --length 20 --q 20 --fastqc -o trim/ sample.fastq.gz

# spingo
########
gunzip -c sample.fq.gz | sed -n '1~4s/^@/>/p;2~4p' > sample.fa
spingo -d ~/spingo/database/RDP_11.2.species.fa -i sample.fa -p 10 > spingo/sample.out

~/src/spingo/spingo_summary sample.out -l 3 -s 0.5 -t 0.7 -p > species_rel/sample.txt
~/src/spingo/spingo_summary sample.out -l 3 -s 0.5 -t 0.7 > species_abs/sample.txt



##############
# paired-end
##############

# QC
####
trim_galore --cores 8 --paired --length 20 --q 20 --fastqc -o trim/ sample_1.fastq.gz sample_2.fastq.gz


# spingo
########
gunzip -c sample_1_val_1.fq.gz | sed -n '1~4s/^@/>/p;2~4p' > sample.1
gunzip -c sample_2_val_2.fq.gz | sed -n '1~4s/^@/>/p;2~4p' > sample.2

cat sample.1 sample.2 > sample.fa
spingo -d ~/spingo/database/RDP_11.2.species.fa -i sample.fa -p 10 > spingo/sample.out

~/src/spingo/spingo_summary sample.out -l 3 -s 0.5 -t 0.7 -p > species_rel/sample.txt
~/src/spingo/spingo_summary sample.out -l 3 -s 0.5 -t 0.7 > species_abs/sample.txt




##########################################
# shotgun metagenome
##########################################

# QC
####
trim_galore --cores 8 --paired --length 20 --q 20 --fastqc -o trim/ sample_1.fastq.gz sample_2.fastq.gz


# filter human and phiX seqs
##############################

for f1 in *_1_val_1.fq.gz
do
f2=${f1%%_1_val_1.fq.gz}"_2_val_2.fq.gz"
name=$(basename $f1 _1_val_1.fq.gz)

#download kraken2_human_and_phiX_db: https://doi.org/10.6084/m9.figshare.12525242.v1
kraken2 --db ~/kraken2_human_and_phiX_db/ --threads 10 --report sample.Host.report --gzip-compressed --unclassified-out sample.noHost#.fq --paired sample_1_val_1.fq.gz sample_2_val_2.fq.gz > sample.kraken.Host.out
pigz -p 10 sample.noHost_1.fq
pigz -p 10 sample.noHost_2.fq

# metaphlan3
##############
metaphlan sample.noHost_1.fq.gz,sample.noHost_2.fq.gz --nproc 10 --bowtie2out sample.bt2.bz2 --input_type fastq -t rel_ab_w_read_stats --bowtie2db ~/humann_db/mpa_v30_CHOCOPhlAn_201901 -o sample.txt

# humann3
###############
cat sample.noHost_1.fq sample.noHost_2.fq > sample.fastq
humann --resume --input sample.fastq --output humann_output/ --threads 10 --nucleotide-database ~/humann_db/chocophlan --protein-database ~/humann_db/uniref --metaphlan-options "--bowtie2db ~/humann_db/mpa_v30_CHOCOPhlAn_201901"

# regroup functional potential
##############################
humann_join_tables --input genefamilies --output sample_genefamilies.tsv --file_name genefamilies
humann_renorm_table --input sample_genefamilies.tsv --output sample_genefamilies_rel.tsv --units relab

humann_regroup_table --input sample_genefamilies.tsv --output sample_ko.tsv --custom ~/humann_db/utility_mapping/map_ko_uniref90.txt.gz
humann_regroup_table --input sample_genefamilies_rel.tsv --output sample_ko_rel.tsv --custom ~/humann_db/utility_mapping/map_ko_uniref90.txt.gz

humann_regroup_table --input sample_genefamilies.tsv --output sample_ec.tsv --custom ~/humann_db/utility_mapping/map_level4ec_uniref90.txt.gz
humann_regroup_table --input sample_genefamilies_rel.tsv --output sample_ec_rel.tsv --custom ~/humann_db/utility_mapping/map_level4ec_uniref90.txt.gz
humann_rename_table --input sample_ec.tsv --output sample_ec_named.tsv --names ec
humann_rename_table --input sample_ec_rel.tsv --output sample_ec_rel_named.tsv --names ec

humann_join_tables --input pathabundance --output sample_pathabundance.tsv --file_name pathabundance
humann_renorm_table --input sample_pathabundance.tsv --output sample_pathabundance_rel.tsv --units relab

# CAZyme mapping file: CAZyme-mapping-database/map_cazy_uniref90.txt
humann_regroup_table --input sample_genefamilies.tsv --output sample_cazy.tsv --custom ~/humann_db/utility_mapping/map_cazy_uniref90.txt.gz
humann_regroup_table --input sample_genefamilies_rel.tsv --output sample_cazy_rel.tsv --custom ~/humann_db/utility_mapping/map_cazy_uniref90.txt.gz

