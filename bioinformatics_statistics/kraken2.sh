#######################################
# Building a custom kraken2 database 
# https://ccb.jhu.edu/software/kraken/MANUAL.html#custom-databasess
#######################################
kraken2-build --download-taxonomy --db kraken_2022_02_18 --threads 30
kraken2-build --download-library bacteria --db kraken_2022_02_18 --threads 30
kraken2-build --download-library fungi --db kraken_2022_02_18 --threads 30
kraken2-build --download-library archaea --db kraken_2022_02_18 --threads 30
kraken2-build --download-library viral --db kraken_2022_02_18 --threads 30
kraken2-build --download-library protozoa --db kraken_2022_02_18 --threads 30


##############################
# run kraken2
##############################
kraken2 --db kraken_2022_02_18 --paired sample.noHost_1.fq sample.noHost_2.fq --threads 30 --use-names --report kraken_report --report-zero-counts --output kraken.out


##############################
# parse report from kraken2
##############################
#https://github.com/npbhavya/Kraken2-output-manipulation/blob/master/kraken-multiple-taxa.py

# Kingdom level
python kraken-multiple-taxa.py -d report/ -r K -c 1 -o report/temp
sed -e "s/\[//g;s/\]//g;s/'//g;s|\t|,|g" report/temp | sed 's@report/@@g' | sed 's/.tax.report//g' > kingdom_rel.csv
rm report/temp


# Species level
python kraken-multiple-taxa.py -d report/ -r S -c 1 -o report/temp
sed -e "s/\[//g;s/\]//g;s/'//g;s|\t|,|g" report/temp | sed 's@report/@@g' | sed 's/.tax.report//g' > species_rel.csv
rm report/temp

