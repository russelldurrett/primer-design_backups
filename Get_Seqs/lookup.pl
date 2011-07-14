#!usr/bin/perl


#script to mine gene CDSs using Gene Names

use Bio::SeqIO;

print "\nWhich Database (E.coli or D.rad)?\n";
my $organism = <STDIN>;
chomp $organism;

my @db_files;

if ($organism =~ /coli/){
	@db_files = ('./genomes/NC_000913.gbk');
	$tag_name = 'gene';
} elsif ($organism =~ /rad/){
	@db_files = ('./genomes/NC_000958.gbk', './genomes/NC_000959.gbk', './genomes/NC_000959.gbk', './genomes/NC_001263.gbk', './genomes/NC_001264.gbk') ;
	$tag_name = 'locus_tag';
} else {print "couldn't understand...exiting\n\n"; exit;}



#Store genes in Hash to reference at every gene name more quickly
my %genes;
open (GENES, 'genes.txt');
while (<GENES>){
	chomp;
	$gene_name = uc $_;
	$genes{$gene_name}++
}
close GENES;

print "searching for these genes:\n";
foreach $gene_key (keys %genes){
	print $gene_key . "\n"
}



print "\n\nparsing @db_files\n.........\n";

open (OUT, '>geneCDSs.fa');

foreach $db_file (@db_files){

	my $seqio_object = Bio::SeqIO->new(-file => $db_file);
	my $seq_object = $seqio_object->next_seq;
	 
	for my $feat_object ($seq_object->get_SeqFeatures) {
		if ($feat_object->primary_tag eq "CDS") { 
	   		if ($feat_object->has_tag($tag_name)) {
	      		#print $feat_object->get_tag_values('gene');
	      		#print "\n";
	      		for $this_gene ($feat_object->get_tag_values($tag_name)) {
	      			$this_gene_name = uc $this_gene;
	      			if (exists $genes{$this_gene_name}) {
	      				print $this_gene . " FOUND\n";
	      				print $feat_object->seq->seq(),"\n";
	      				my $bp = length $feat_object->seq->seq();
	      				print OUT ">$this_gene $bp" . "bp\n";
	      				print OUT $feat_object->seq->seq() . "\n";
	      				delete($genes{$this_gene_name});
	     	  		}
	     		}
	      	}
		}
	}
}


print "\nThese weren't found in the genome file:\n";
foreach $gene_key (keys %genes){
	print $gene_key . "\n"
}