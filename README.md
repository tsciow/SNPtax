# SNPtax
## A pipeline to identify SNPs specific for selected taxa

This pipeline consists of 2 perl-scripts which require bioperl.
Prank (http://wasabiapp.org/software/prank/) is used to align sequences.

### Requires:
- BioPerl
- Getopt::Std
- Try::Tiny


Please keep in mind that these tools **DO NOT** use the NCBI taxonomy database or similar, but rather use a quick and dirty hack.
The last common ancestor is determined simply by string comparison of the taxonomy and taking the longest common prefix of 2 strings.
For example :
- Eukaryota_Viridiplantae_Streptophyta_Embryophyta_Tracheophyta_Spermatophyta_Magnoliophyta_eudicotyledons_Gunneridae_

and
- Eukaryota_Viridiplantae_Streptophyta_Embryophyta_Tracheophyta_Spermatophyta_Pinidae_Cupressales_Taxaceae_Taxus_Taxus_baccata_

will result in 

- Eukaryota_Viridiplantae_Streptophyta_Embryophyta_Tracheophyta_Spermatophyta_

as the common taxonomy.

This can lead to some strange artefacts, for example
- eudicotyledons,Gunneridae,Pentapetalae,rosids,malvids,Malvales,Malvaceae,Bombacoideae,Bombax,Bombax_ceiba,
- eudicotyledons,Gunneridae,Pentapetalae,rosids,malvids,Myrtales,Myrtaceae,Myrtoideae,Eucalypteae,Eucalyptus,Eucalyptus_grandis,

will result in 
- eudicotyledons,Gunneridae,Pentapetalae,rosids,malvids,Myrtales,M

obviously not the real last common ancestor but the common prefix of the two taxonomy strings

Just keep this in mind when browsing/screening the result files for the taxa of your interest!


Another thing to consider:
for our pipeline we used Prank as alignment tool.
In fasta headers Prank replaces commata "," with underscores "_".
Other alignment tools might behave differently, so you'd have to adjust the screening step later on.

The positions in the alignments are given as absolute positions and as relative postitions realtive to two refernece-species.
These reference species are hard-coded in the script and need to be changed there. Don't forget the terminal delimiter (in our case with prank aligner the terminal underscore _ , for example "Helicobacter_pylorii_").



### SNPtax_extract_genes_from_gbk.pl
- INPUT: a (organellar) genome in Genbank format
- OUTPUT: several single-entry fasta files, one for each annotated gene from the Genbank-input (GeneID followed by a double underscore __ followed by the complete taxonomy separated by commata ,
- OUTPUT to STDOUT: a list of all genes annotated in the Genbank-input


### SNPtax_process_alignment.pl
- INPUT: a multiple alignment file in fasta format
- OUTPUT: 3 files:
  - *.used_taxa: a file listing all taxa present in the alignemnt file
  - *.SNP2: contains, for each heterogenous position (absolute) in the alignment, the bases occuring at that position and the "last common ancestor" carrying given base
  - *.SNP: The most interesting output file: lists all positions with a taxon-specific base together with the taxon this base is specific for.
                Screening this file for our taxa of interest will give us a list of taxon-specific SNPs
                

## Example
Let's assume we have four chloroplast-genomes:
- NC_027425.1.gbk  Populus tremula
- NC_036929.1.gbk  Fagus engleriana
- NC_041252.1.gbk  Fagus crenata
- NC_041437.1.gbk  Fagus sylvatica

We then run SNPtax_extract_genes_from_gbk.pl on each of these Genbank files:
- /path/to/SNPtax_extract_genes_from_gbk.pl -i NC_027425.1.gbk
- /path/to/SNPtax_extract_genes_from_gbk.pl -i NC_036929.1.gbk
- /path/to/SNPtax_extract_genes_from_gbk.pl -i NC_041252.1.gbk
- /path/to/SNPtax_extract_genes_from_gbk.pl -i NC_041437.1.gbk

This will give us (among many other files) 4 fasta files for the psac-Gene from these 4 genomes
- psac__Eukaryota,Viridiplantae,Streptophyta,Embryophyta,Tracheophyta,Spermatophyta,Magnoliophyta,eudicotyledons,Gunneridae,Pentapetalae,rosids,fabids,Fagales,Fagaceae,Fagus,Fagus_crenata,.fna
- psac__Eukaryota,Viridiplantae,Streptophyta,Embryophyta,Tracheophyta,Spermatophyta,Magnoliophyta,eudicotyledons,Gunneridae,Pentapetalae,rosids,fabids,Fagales,Fagaceae,Fagus,Fagus_engleriana,.fna
- psac__Eukaryota,Viridiplantae,Streptophyta,Embryophyta,Tracheophyta,Spermatophyta,Magnoliophyta,eudicotyledons,Gunneridae,Pentapetalae,rosids,fabids,Fagales,Fagaceae,Fagus,Fagus_sylvatica,.fna
- psac__Eukaryota,Viridiplantae,Streptophyta,Embryophyta,Tracheophyta,Spermatophyta,Magnoliophyta,eudicotyledons,Gunneridae,Pentapetalae,rosids,fabids,Malpighiales,Salicaceae,Saliceae,Populus,Populus_tremula,.fna

We then concatenate all psac-sequences into one file:
- cat psac*.fna > psac.fasta

Then we run the prank aligner on this multi-fasta file:
- /path/to/prank -d=psac.fasta -o=psac.aligned.fasta 

Finally we run SNPtax_process_alignment.pl on the aligned fasta file (in this example we have 4 species in the alignment, so let's just set n to a minimum of 2 species; this parameter is more useful for the automated analysis of larger numbers of taxa)
- /path/to/SNPtax_process_alignment.pl -i psac.aligned.fasta -n 2

This will result in 3 files:
- psac.aligned.fasta.best.fas.SNP
- psac.aligned.fasta.best.fas.SNP2
- psac.aligned.fasta.best.fas.used_taxa

We'll have a closer look at psac.aligned.fasta.best.fas.SNP
```
Abs.:18	Populus tremula:18	Pinus taeda:1	Spermatophyta_Magnoliophyta_eudicotyledons_Gunneridae_Pentapetalae_rosids_fabids_Malpighiales_Salicaceae_Saliceae_Populus_Populus_tremula_		A
Abs.:18	Populus tremula:18	Pinus taeda:1	Spermatophyta_Magnoliophyta_eudicotyledons_Gunneridae_Pentapetalae_rosids_fabids_Fagales_Fagaceae_Fagus_Fagus_		G
Abs.:21	Populus tremula:21	Pinus taeda:1	Spermatophyta_Magnoliophyta_eudicotyledons_Gunneridae_Pentapetalae_rosids_fabids_Malpighiales_Salicaceae_Saliceae_Populus_Populus_tremula_		C
Abs.:21	Populus tremula:21	Pinus taeda:1	Spermatophyta_Magnoliophyta_eudicotyledons_Gunneridae_Pentapetalae_rosids_fabids_Fagales_Fagaceae_Fagus_Fagus_		T
Abs.:55	Populus tremula:55	Pinus taeda:1	Spermatophyta_Magnoliophyta_eudicotyledons_Gunneridae_Pentapetalae_rosids_fabids_Malpighiales_Salicaceae_Saliceae_Populus_Populus_tremula_		C
Abs.:55	Populus tremula:55	Pinus taeda:1	Spermatophyta_Magnoliophyta_eudicotyledons_Gunneridae_Pentapetalae_rosids_fabids_Fagales_Fagaceae_Fagus_Fagus_		A
Abs.:57	Populus tremula:57	Pinus taeda:1	Spermatophyta_Magnoliophyta_eudicotyledons_Gunneridae_Pentapetalae_rosids_fabids_Malpighiales_Salicaceae_Saliceae_Populus_Populus_tremula_		G
Abs.:57	Populus tremula:57	Pinus taeda:1	Spermatophyta_Magnoliophyta_eudicotyledons_Gunneridae_Pentapetalae_rosids_fabids_Fagales_Fagaceae_Fagus_Fagus_		A
Abs.:63	Populus tremula:63	Pinus taeda:1	Spermatophyta_Magnoliophyta_eudicotyledons_Gunneridae_Pentapetalae_rosids_fabids_Malpighiales_Salicaceae_Saliceae_Populus_Populus_tremula_		C
Abs.:63	Populus tremula:63	Pinus taeda:1	Spermatophyta_Magnoliophyta_eudicotyledons_Gunneridae_Pentapetalae_rosids_fabids_Fagales_Fagaceae_Fagus_Fagus_		T
Abs.:69	Populus tremula:69	Pinus taeda:1	Spermatophyta_Magnoliophyta_eudicotyledons_Gunneridae_Pentapetalae_rosids_fabids_Malpighiales_Salicaceae_Saliceae_Populus_Populus_tremula_		G
Abs.:69	Populus tremula:69	Pinus taeda:1	Spermatophyta_Magnoliophyta_eudicotyledons_Gunneridae_Pentapetalae_rosids_fabids_Fagales_Fagaceae_Fagus_Fagus_		A
Abs.:99	Populus tremula:99	Pinus taeda:1	Spermatophyta_Magnoliophyta_eudicotyledons_Gunneridae_Pentapetalae_rosids_fabids_Malpighiales_Salicaceae_Saliceae_Populus_Populus_tremula_		G
Abs.:99	Populus tremula:99	Pinus taeda:1	Spermatophyta_Magnoliophyta_eudicotyledons_Gunneridae_Pentapetalae_rosids_fabids_Fagales_Fagaceae_Fagus_Fagus_		A
Abs.:135	Populus tremula:135	Pinus taeda:1	Spermatophyta_Magnoliophyta_eudicotyledons_Gunneridae_Pentapetalae_rosids_fabids_Fagales_Fagaceae_Fagus_Fagus_		G
Abs.:135	Populus tremula:135	Pinus taeda:1	Spermatophyta_Magnoliophyta_eudicotyledons_Gunneridae_Pentapetalae_rosids_fabids_Malpighiales_Salicaceae_Saliceae_Populus_Populus_tremula_		A
Abs.:138	Populus tremula:138	Pinus taeda:1	Spermatophyta_Magnoliophyta_eudicotyledons_Gunneridae_Pentapetalae_rosids_fabids_Fagales_Fagaceae_Fagus_Fagus_		G
Abs.:138	Populus tremula:138	Pinus taeda:1	Spermatophyta_Magnoliophyta_eudicotyledons_Gunneridae_Pentapetalae_rosids_fabids_Malpighiales_Salicaceae_Saliceae_Populus_Populus_tremula_		A
Abs.:150	Populus tremula:150	Pinus taeda:1	Spermatophyta_Magnoliophyta_eudicotyledons_Gunneridae_Pentapetalae_rosids_fabids_Fagales_Fagaceae_Fagus_Fagus_		T
Abs.:150	Populus tremula:150	Pinus taeda:1	Spermatophyta_Magnoliophyta_eudicotyledons_Gunneridae_Pentapetalae_rosids_fabids_Malpighiales_Salicaceae_Saliceae_Populus_Populus_tremula_		C
Abs.:156	Populus tremula:156	Pinus taeda:1	Spermatophyta_Magnoliophyta_eudicotyledons_Gunneridae_Pentapetalae_rosids_fabids_Fagales_Fagaceae_Fagus_Fagus_sylvatica_		A
Abs.:162	Populus tremula:162	Pinus taeda:1	Spermatophyta_Magnoliophyta_eudicotyledons_Gunneridae_Pentapetalae_rosids_fabids_Malpighiales_Salicaceae_Saliceae_Populus_Populus_tremula_		C
Abs.:162	Populus tremula:162	Pinus taeda:1	Spermatophyta_Magnoliophyta_eudicotyledons_Gunneridae_Pentapetalae_rosids_fabids_Fagales_Fagaceae_Fagus_Fagus_		T
Abs.:180	Populus tremula:180	Pinus taeda:1	Spermatophyta_Magnoliophyta_eudicotyledons_Gunneridae_Pentapetalae_rosids_fabids_Malpighiales_Salicaceae_Saliceae_Populus_Populus_tremula_		T
Abs.:180	Populus tremula:180	Pinus taeda:1	Spermatophyta_Magnoliophyta_eudicotyledons_Gunneridae_Pentapetalae_rosids_fabids_Fagales_Fagaceae_Fagus_Fagus_		G
Abs.:192	Populus tremula:192	Pinus taeda:1	Spermatophyta_Magnoliophyta_eudicotyledons_Gunneridae_Pentapetalae_rosids_fabids_Fagales_Fagaceae_Fagus_Fagus_		C
Abs.:192	Populus tremula:192	Pinus taeda:1	Spermatophyta_Magnoliophyta_eudicotyledons_Gunneridae_Pentapetalae_rosids_fabids_Malpighiales_Salicaceae_Saliceae_Populus_Populus_tremula_		T
Abs.:205	Populus tremula:205	Pinus taeda:1	Spermatophyta_Magnoliophyta_eudicotyledons_Gunneridae_Pentapetalae_rosids_fabids_Fagales_Fagaceae_Fagus_Fagus_		T
Abs.:205	Populus tremula:205	Pinus taeda:1	Spermatophyta_Magnoliophyta_eudicotyledons_Gunneridae_Pentapetalae_rosids_fabids_Malpighiales_Salicaceae_Saliceae_Populus_Populus_tremula_		C
Abs.:234	Populus tremula:234	Pinus taeda:1	Spermatophyta_Magnoliophyta_eudicotyledons_Gunneridae_Pentapetalae_rosids_fabids_Fagales_Fagaceae_Fagus_Fagus_		C
Abs.:234	Populus tremula:234	Pinus taeda:1	Spermatophyta_Magnoliophyta_eudicotyledons_Gunneridae_Pentapetalae_rosids_fabids_Malpighiales_Salicaceae_Saliceae_Populus_Populus_tremula_		T
Abs.:237	Populus tremula:237	Pinus taeda:1	Spermatophyta_Magnoliophyta_eudicotyledons_Gunneridae_Pentapetalae_rosids_fabids_Malpighiales_Salicaceae_Saliceae_Populus_Populus_tremula_		A
Abs.:237	Populus tremula:237	Pinus taeda:1	Spermatophyta_Magnoliophyta_eudicotyledons_Gunneridae_Pentapetalae_rosids_fabids_Fagales_Fagaceae_Fagus_Fagus_		C
Abs.:245	Populus tremula:245	Pinus taeda:1	Spermatophyta_Magnoliophyta_eudicotyledons_Gunneridae_Pentapetalae_rosids_fabids_Malpighiales_Salicaceae_Saliceae_Populus_Populus_tremula_		A
Abs.:245	Populus tremula:245	Pinus taeda:1	Spermatophyta_Magnoliophyta_eudicotyledons_Gunneridae_Pentapetalae_rosids_fabids_Fagales_Fagaceae_Fagus_Fagus_		G
```
##### Note the tailing underscore followed by a TAB followed by the Nucleotide ! We'll use this to screen the result file
  
By using the grep command we can screen psac.aligned.fasta.best.fas.SNP for the taxa of our interest:
- For SNPs common and specific to all Fagus species we'd use (Note the _\t )
```
grep -P "Fagus_\t" psac.aligned.fasta.best.fas.SNP
```
- For SNPs specific for Fagus sylvatica we'd use (Note the _\t )
```
grep -P "Fagus_sylvatica_\t" psac.aligned.fasta.best.fas.SNP
```

##### Note:
#####	prank aligner replaces , with _
#####	the taxonomy-string always ends with a _
#####	in the SNP files the columns are delimited by a TAB ,
#####	so we'll search the SNP files for TaxonOfInterest_TAB , for example "Fagaceae_    "
#####	we'll need to use grep with the -P (use perl syntax, so it interprets \t as TAB)






## Pipeline:
- 1.) run SNPtax_extract_genes_from_gbk.pl on the Genbank files of the taxa of interest to extract genic sequences
        - a list of all gene names will be written to STDOUT
        - the different taxpnomic levels taken from the ORGANISM flag in the Genbank file will be comma-delimited; 
          - also at the end a comma will be inserted
          - prank aligner will replace thes commata with an underscore character _
- 2.) for each gene concatenate all sequence files into one file
- 3.) run prank to align the gene sequences
- 4.) run SNPtax_process_alignment.pl to process the alignment files
          - parameter -n is the minimum number of entrie in the alignment needed so the program will continue the analysis
            (if you' re looking for SNPs/Markers to differentiate 50 different species, there's no point in looking
              at genes only present in a few species)
- 5.) browse/screen the output files for the taxa of your interest
        - the SNP files will be the most interesting files
          - when screening these for taxa of interest, don't forget the terminal underscore character _
            for example: grep -P "Fagaceae_\t" atpA__.prankalignment.SNP
            
  
  

A batch script could look something like this:

```
#!/bin/bash

#1.) Extract gene sequences from all Genbank files in this folder

for p in ./1_gbk/*.gbk; do ./SNPtax_extract_genes_from_gbk.pl -i $p; done | sort -u > genes.lst

wait

#the names of all genes present in ANY of the GBK files is written to genes.lst


#2.) move all fna files to one temporary directory
mkdir 2_genes_fasta
mv *.fna ./2_genes_fasta/

#3.) for each gene present in any of the genbank files concatenate all fna files into one multi-fasta file


while read l; do cat ./2_genes_fasta/$l*.fna > $l.fasta ; done < genes.lst

wait
mv *.fasta ./3_concatenated_fasta/

#3. for each gene (multi-fasta file) run prank aligner

while read f;do /change/path/bin/prank -d=./3_concatenated_fasta/$f.fasta -o=./4_prank_alignments/$f ;done < genes.lst

wait

#4. for each alignment run our SNP-calling script

while read a; do ./SNPtax_process_alignment.pl -i $a*.fas -n 10 ; done < genes.lst

wait

#6. Parse the SNP files for taxon specific SNPs

#Note:
#	prank aligner replaces , with _
#	the taxonomy-string always ends with a _
#	in the SNP files the columns are delimited by a <TAB> ,
#	so we'll search the SNP files for TaxonOfInterest_<TAB> , for example "Fagaceae_	"
#	we'll need to use grep with the -P (use perl syntax, so it interprets \t as <TAB>)

grep -P "Pinaceae_\t" *.SNP > Marker_Pinaceae.txt

#or

grep -P "Fagus_sylvatica_\t" *.SNP > Marker_Fagus_sylvatica.txt



