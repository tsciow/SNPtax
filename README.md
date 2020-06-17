# SNPtax
A pipeline to identify SNPs specific for selected taxa

This pipeline consists 2 perl-scripts which require bioperl.
Prank (http://wasabiapp.org/software/prank/) is used to align sequences.

Requires:
BioPerl
Getopt::Std
Try::Tiny


Please keep in mind that these tools DO NOT use the NCBI taxonomy database or similar, but rather uses a quick and dirty hack.
The last common ancestor is determined simply by string comparison of the taxonomy and taking the longest common prefix of 2 strings.
For example :
Eukaryota_Viridiplantae_Streptophyta_Embryophyta_Tracheophyta_Spermatophyta_Magnoliophyta_eudicotyledons_Gunneridae_
and
Eukaryota_Viridiplantae_Streptophyta_Embryophyta_Tracheophyta_Spermatophyta_Pinidae_Cupressales_Taxaceae_Taxus_Taxus_baccata_

will result in 
Eukaryota_Viridiplantae_Streptophyta_Embryophyta_Tracheophyta_Spermatophyta_
as the common taxonomy.

This can lead to some strange artefacts, for example
>eudicotyledons,Gunneridae,Pentapetalae,rosids,malvids,Malvales,Malvaceae,Bombacoideae,Bombax,Bombax_ceiba,
>eudicotyledons,Gunneridae,Pentapetalae,rosids,malvids,Myrtales,Myrtaceae,Myrtoideae,Eucalypteae,Eucalyptus,Eucalyptus_grandis,
will result in 
>eudicotyledons,Gunneridae,Pentapetalae,rosids,malvids,Myrtales,M
obviously not the real last common ancestor but the common prefix of the two taxonomy strings

Just keep this in mind when browsing/screening the result files for the taxa of your interest!


Another thing to consider:
for our pipeline we used Prank as alignment tool.
Prank replaces commata "," with underscores "_".
Other alignment tools might behave differently, so you'd have to adjust the screening step later on.



