#! /usr/bin/perl
use strict;
use Bio::SeqIO;
use Getopt::Std;
use Bio::SimpleAlign;
use Bio::AlignIO;

sub usage () {
	print << "EOF";
	SNPtax_process_alignment.pl:
	This Program reads a fasta-alignment file and detects, for each position in the alignment, the taxa for which a nucleotide is characteristic.
	Three output files will be generated:
	1.) *.used_taxa : a list of all taxa present in the alignment file
	2.) *.SNP2 : for each position in the alignment the occuring bases with the corresponding last common ancestor
	3.) *.SNP :  
	Parameters:
	-i: input fasta-alignment file
	-n: minimum number of entries/taxa/species required to continue analysis
	Example:
	SNPtax_process_alignment.pl -i petD.fasta -n 10
  
EOF
	exit;
};

sub longest_common_prefix {# truncated up to the last comma
    # longest_common_prefix( $|@ ): returns $
    # URLref: http://linux.seindal.dk/2005/09/09/longest-common-prefix-in-perl
    # find longest common prefix of scalar list
    my $prefix = shift;
    for (@_) {
        chop $prefix while (! /^\Q$prefix\E/);
        }
    if ((substr $prefix, -1) ne ","){
       my $i = rindex($prefix,"_");
       $prefix = (substr ($prefix,0,$i))."_";
    }
    return $prefix;
}



if ($#ARGV < 1){
	usage();
}

my %opts;
getopts ('i:n:', \%opts) or usage();

my $infile = $opts{i};
my $min_number = $opts{n};
my $ali_in = Bio::AlignIO->new('-file' => $infile , '-format' => 'fasta');

my $organism;
my $taxonomy;


open (SNP,'>',$infile.".SNP"); # List of Golden SNPs: Absolute Position in alignment; Position in Polpulus tremula gene; Poition in Pinus taeda gene; Tax. level for which this Posiiton carries a golden SNP; nucleotide spicific for this taxon
open (SNP2,'>',$infile.".SNP2"); # List of nucleotides present in each position with the taxonomy for which this base is specific
open (TAX,'>',$infile.".used_taxa");# List of taxa present in the alignment (sometimes a gene is not annotated in the published genbank file, so some species/gene combinations might be missing 

my $aln = $ali_in->next_aln();

if ($aln->num_sequences < $min_number){#if less than $min_number ( -n parameter) Sequences in Alignment, quit. No point in picking SNPs 
            
close (SNP);
close (SNP2);
close (TAX);
      
die;
}



foreach my $point (0..$aln->length()-1){
my %base;#hash base->longest common substring of taxonomy
my %taxon;
my $Laub_position=0;
my $Nadel_position=0;

my $gaps=0;#count the number of gaps in a position
my $nuc=0;# number of nucleotides in a position of the alignment; used to check if only one seuqence actually has a base and all others have gap in this position. If so, ignore it

                                                                                                             
                                                                                                                
foreach my $seq ( $aln->each_seq() ) {
if (substr($seq->seq,$point,1) eq "-"){
   $gaps++;
}
else {
     $nuc++;
}




#Positions relative to reference genomes        
        if (index($seq->id(),"Populus_tremula_") != -1){
                if (defined $seq->location_from_column($point+1)){
                   $Laub_position = ($seq->location_from_column($point+1)->start)-1;
                }
        }
        if (index($seq->id(),"Pinus_taeda_") != -1){
           if (defined $seq->location_from_column($point+1)){
              $Nadel_position = ($seq->location_from_column($point+1)->start)-1;
           } 
        }
        
#fill hash %base: key= nucleotide at that posiiton, value=last common ancestor that has this base in this position
        if (defined $base{substr($seq->seq,$point,1)}){
           my @tmparray = ($base{substr($seq->seq,$point,1)}, substr($seq->id(),62));
           $base{substr($seq->seq,$point,1)} = longest_common_prefix (@tmparray);
        }
        else{
             $base{substr($seq->seq,$point,1)} = substr($seq->id(),62);
        }
        
#fill hash %taxon: key: allcombinations of the given taxonomic levels (i.e. Embryophyta,Tracheophyta,Spermatophyta: this results in the keys Embryophyta; Embryophyta,Tracheophyta and Embryophyta,Tracheophyta,Spermatophyta) ; value = the nucleotide in this position. if it differs from the base stored in the hash it'll be set to "Not Golden";

        my @felder = split /_/, substr($seq->id(),62);# clip of common part of taxonomy: remove Eukaryota,Viridiplantae,Streptophyta,
        my $anzahl=@felder;        
        my $tax="";
        for my $i ( 0 .. $anzahl-1){# -1 because the last letter is ','
            $tax.=$felder[$i]."_"; # add an underscore _ as last letter; we'll need this later when processing the SNP-files
            if (defined $taxon{$tax}){
               if ($taxon{$tax} ne substr($seq->seq,$point,1)){
                  $taxon{$tax}="Not Golden";
               }            
            }
            elsif (length($tax)>5) {
                 $taxon{$tax}=substr($seq->seq,$point,1);
            }
        }                                                                                              
}

my $key;
my $value;
if (keys %base > 1){
print SNP2 "\n------------------------------------------------------------\n",$point+1;
while (($key, $value) = each (%base)){
      
         if ($value ne ''){
            print SNP2 "\n$value\t\t$key";
         }        
      }
}


if (keys %base > 1){
#print SNP "\n------------------------------------------------------------\n$point";
while (($key, $value) = each (%base)){
      if($nuc > 1){ # filter out those SNPs, where only one Sequence has a base and all others have a gap
         if ($value ne ''){
            if (($taxon{$value} ne 'Not Golden')&&(length($value)>2)){ #&&(length($value)>2 not necessary anymore; due to tax. misannotations in genbank files we had some artifacts here for the common taxonomy $value (e.g. Mag_ because of Magnoliophyta and Magnoliopsida
               print SNP "Abs.:",$point+1,"\tPopulus tremula:",$Laub_position+1,"\tPinus taeda:",$Nadel_position+1,"\t$value\t\t$key\n";
            }
         }        
      }
      }
}
        
}

##################################################################################
#print list of species that where taken into account/present in alignment

print TAX "Sequences taken into account:\n\n";
my $seq;
foreach $seq ($aln->each_seq()){
      print TAX "\n".$seq->id();
}


            
close (SNP);
close (SNP2);
close (TAX);
      
      
