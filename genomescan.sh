#!/bin/bash 

#This shell script performs some Fst/popgen outler analysis on a vcf file. 
#it was inspired by speciationgenomics.github.io/per_site_Fst/ and the sliding window analysis

#extract the population on which you want to calculate the pers ite fst 

~/Programs/bin/bcftools query -l ~/SouthCenter.snps.vcf.gz | grep "SC_4-8" > four

#repeat for the the other popualtions:

~/Programs/bin/bcftools query -l ~/SouthCenter.snps.vcf.gz| grep "SC_lake" > lake
~/Programs/bin/bcftools query -l ~/SouthCenter.snps.vcf.gz| grep "SC_12-16" > twelve
~/Programs/bin/bcftools query -l ~/SouthCenter.snps.vcf.gz| grep "SC_20-24" > twenty
~/Programs/bin/bcftools query -l ~/SouthCenter.snps.vcf.gz| grep "SC_52-56" > fiftytwo
~/Programs/bin/bcftools query -l ~/SouthCenter.snps.vcf.gz| grep "SC_60-64" > sixty

# Do the calculations using vcf tools 

~/Programs/bin/vcftools --gzvcf ~/SouthCenter.snps.vcf.gz \
--weir-fst-pop four \
--weir-fst-pop twelve \
--weir-fst-pop twenty \
--weir-fst-pop fiftytwo \
--weir-fst-pop sixty \
--out SouthCenter
# Use r to parse the output see the Popgenomics.r script for details 

#Generate the sliding window stats 
# Convert to the .geno formats

python2 ~/Programs/genomics_general/VCF_processing/parseVCF.py -i ~/SouthCenter.snps.vcf.gz -o SouthCenter.geno.gz

# Generate the pops file 

~/Programs/bin/bcftools query -l ~/SouthCenter.snps.vcf.gz > individuals.txt

cat individuals.txt |while read line; do echo $line| tr "\n" "\t"; echo $line|awk -F "_" '{print $2}'; done >pops.file

#Calc the sliding window stats 
python2 ~/Programs/genomics_general/popgenWindows.py -g SouthCenter.geno.gz -o SouthCenter.window.csv.gz -f phased -w 10000 -m 3 -s 2000 -p 4-8 -p 12-16 -p 20-24 -p lake -p 52-56 -p 60-64 --popsFile pops.file

#parse file in R with popgenomics.r 

#done

