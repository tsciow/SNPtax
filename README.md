# SNPtax
##A pipeline to identify SNPs specific for selected taxa

This pipeline consists 2 perl-scripts which require bioperl.
Prank (http://wasabiapp.org/software/prank/) is used to align sequences.

###Requires:
BioPerl
Getopt::Std
Try::Tiny


Please keep in mind that these tools **DO NOT** use the NCBI taxonomy database or similar, but rather uses a quick and dirty hack.
The last common ancestor is determined simply by string comparison of the taxonomy and taking the longest common prefix of 2 strings.
For example :
Eukaryota_Viridiplantae_Streptophyta_Embryophyta_Tracheophyta_Spermatophyta_Magnoliophyta_eudicotyledons_Gunneridae_
and
Eukaryota_Viridiplantae_Streptophyta_Embryophyta_Tracheophyta_Spermatophyta_Pinidae_Cupressales_Taxaceae_Taxus_Taxus_baccata_

will result in 
Eukaryota_Viridiplantae_Streptophyta_Embryophyta_Tracheophyta_Spermatophyta_
as the common taxonomy.

This can lead to some strange artefacts, for example
eudicotyledons,Gunneridae,Pentapetalae,rosids,malvids,Malvales,Malvaceae,Bombacoideae,Bombax,Bombax_ceiba,
eudicotyledons,Gunneridae,Pentapetalae,rosids,malvids,Myrtales,Myrtaceae,Myrtoideae,Eucalypteae,Eucalyptus,Eucalyptus_grandis,
will result in 
eudicotyledons,Gunneridae,Pentapetalae,rosids,malvids,Myrtales,M
obviously not the real last common ancestor but the common prefix of the two taxonomy strings

Just keep this in mind when browsing/screening the result files for the taxa of your interest!


Another thing to consider:
for our pipeline we used Prank as alignment tool.
Prank replaces commata "," with underscores "_".
Other alignment tools might behave differently, so you'd have to adjust the screening step later on.




Pipeline:
1.) run SNPtax_extract_genes_from_gbk.pl on the Genbank files of the taxa of interest to extract genic sequences
2.) for each gene concatenate all sequence files into one file
3.) run prank to align the gene sequences
4.) run SNPtax_process_alignment.pl to process the alignment files
5.) browse/screen the output files for the taxa of your interest



#!/bin/bash

#1.) Extract gene sequences from all Genbank files in this folder
mkdir temp_fasta

for p in *.gbk; do ../vTI_extract_gene_nt_seqs_from_gbk.pl -i $p; done | sort -u > genes.lst

wait

#2.) for each gene present in any of the genbank files concatenate all fna files into one multi-fasta file

while read l;
do
cat ./temp_fasta/$l*.fna > $l.fasta
done < genes.lst

wait

#3. for each gene (multi-fasta file) run prank aligner

while read f;
do
/home/schott/bin/prank/bin/prank -d=$f.fasta -o=$f
done < genes.lst

wait

#4. for each alignment run our SNP-calling script

while read a;
do
../vTI_get_SNPS_by_taxonomy_from_alignment.plg -i $a*.fas
done < genes.lst

wait

#5. clean up, delete gene_fasta file
rm -Rf ./temp_fasta/
rm *.fasta



