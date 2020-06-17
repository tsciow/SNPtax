#! /usr/bin/perl
use strict;
use Bio::SeqIO;
use Getopt::Std;
use Try::Tiny;

sub usage () {
	print << "EOF";
	process_gbk.pl:
	this program converts a (multi-)-genbank file.
	The unspliced (can be changed, please see comment) gene sequences of all CDSs are written to separate fasta-files with the complete 
	taxonomy (as it appears in the genbank file under the ORGANSIM tag) as fasta header. The file name consist of the gene-name and the 
	complete taxonomy with a double underscore __ as delimiter.

	Parameters:
	-i: input Genbank-file
	Example:
	SNPtax_extract_genes_from_gbk.pl -i Pinus_taeda.gbk
	reads Pinus_taeda.gbk and writes individual fasta files of all annotated CDS, for example
	sdh3__Eukaryota,Viridiplantae,Streptophyta,Embryophyta,Tracheophyta,Spermatophyta,Pinidae,Pinales,Pinaceae,Pinus,Pinus,Pinus_taeda,.fna

EOF
	exit;
};

if ($#ARGV < 1){
	usage();
}

my %opts;
getopts ('i:', \%opts) or usage();

my $seq_in = Bio::SeqIO->new('-file' => $opts{i}, '-format' => 'genbank');
my $taxonomy;

while ( my $inseq = $seq_in->next_seq ) {
      my @classification = $inseq->species->classification;
      
#@classification contains taxonomy in reverse order (species; genus; family ...), so we need to reverse it      
      my $t;
      foreach $t (reverse @classification){
              $taxonomy = $taxonomy.$t.",";#print "$t;";
              }
      $taxonomy =~ s/ /_/g;# replace " " with "_", so the Species name doesn't contain any white spaces           
	
	for my $feat_object ($inseq->get_SeqFeatures) {
      if ($feat_object->primary_tag eq "CDS") {
      			if ($feat_object->has_tag('gene')) {
         			for my $val ($feat_object->get_tag_values('gene')){
      					my $gene = lc $val =~ s/\//_/r;#make it lower case, so Nad7 and nad7 are equal
					my $seq_out_fna = Bio::SeqIO->new('-file' => ">./".$gene."__".$taxonomy.".fna" , '-format' => 'fasta');
 					try{
           my $seq_obj_fna = Bio::Seq->new(-seq => $feat_object->seq->seq,-display_id => $taxonomy);#In case you want to work with the spliced sequences : $feat_object->spliced_seq->seq,-display_id => $taxonomy);
					$seq_out_fna->write_seq($seq_obj_fna);
     }
     catch{
           warn "Bioperl error";
           };
					print "\n",$gene,"__";
         			}
    			}
   		}
	} 
}




