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



