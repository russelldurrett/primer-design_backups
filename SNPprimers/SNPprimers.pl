#!usr/bin/perl

#SNP Primers Script



# Uses Bioperl to request GenBank files and store in queryable objects. (STREAM_OBJ)

# QUERY IS BASED ON ORGANISM AND GENE FIELDS, BUT CAN EASILYY BE MADE TO SUPPORT OTHERS
use lib '/Volumes/Macintosh HD/Users/Rover/Perl/Master-Primer-Designer/modules';
use lib '/Volumes/Macintosh HD/Users/Rover/Perl/Modules/bioperl';

use Bio::Seq;
use Bio::SeqIO;
use Bio::DB::GenBank;
use Bio::DB::Query::GenBank;
use Extract;

$command_line_input = shift @ARGV;


if ($command_line_input =~ /snp/ && defined($ARGV[0])){
	$testsnp_file = shift @ARGV;
	&snp_primers;
} elsif ($command_line_input =~ /.txt/){
	$testsnp_file = $command_line_input;
	&snp_primers;
} elsif ($command_line_input =~ /A|G|C|T/i){
	$sequence = $command_line_input;
	&command_line_cloning_primers;
} elsif ($command_line_input == ""){
	&genbank_primers;
} else { &print_input_error; }



sub snp_primers{
	
	
	# Written by Russell Durrett in the lab of Chris Mason. rud2004@med.cornell.edu
	
	# BEWARE - because SNPs are sorted by chromosome, the order of SNP primers output is not the same as the input
	
	
	#ASK WHAT GENOME YOU WANT TO WORK WITH
	#DOWNLOAD GENOME DATA IF NOT AVAILABLE
	
	print "What genome are you finding primers for? [Well supported: HG18, HG19 & MM9] \n Genome:  ";
	our $genome = <STDIN>;
	chomp $genome;
	uc $genome;
	
	if ($genome =~ m/^HG18/i){
		#DOWNLOAD GENOME IF HG18 IS SPECIFIED	
		$filename = "./HG18/chr22.fa";
		unless (-e $filename) { 
			print "\n Genome to isolate SNP loci not available, do you want to download Human Genome build 18?   ";
			my $input = "";
			$input = <STDIN>;
			if ($input =~ m/Y/i){
				system "wget --timestamping -P ./HG18/ \'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg18/snp130Mask/*' ";
				system "gunzip ./HG18/*.gz";
				#CURRENTLY NOT SUPPORTING BOWTIE
				#system "unzip *.zip -d bowtie-0.12.7/indexes/";
			} else {print "\n Well, then I can't really do my job. Goodbye. \n\n"; exit; }
		}
	} elsif ($genome =~ m/^HG19/i){
		#DOWNLOAD GENOME IF HG19 IS SPECIFIED	
		$filename = "./HG19/chr22.subst.fa";
		unless (-e $filename) { 
			print "\n Genome to isolate SNP loci not available, do you want to download Human Genome build 19?   ";
			my $input = "";
			$input = <STDIN>;
			if ($input =~ m/Y/i){
				system "wget --timestamping -P ./HG19/ \'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/snp132Mask/*' ";
				system "gunzip ./HG19/*.gz";
			} else {print "\n Well, then I can't really do my job. Goodbye. \n\n"; exit; }
		}
	} elsif ($genome =~ m/^MM9/i){
		#DOWNLOAD GENOME IF MM9 IS SPECIFIED	
		$filename = "./MM9/chr22.subst.fa";
		unless (-e $filename) { 
			print "\n Genome to isolate SNP loci not available, do you want to download Mouse Genome build 9?   ";
			my $input = "";
			$input = <STDIN>;
			if ($input =~ m/Y/i){
				system "wget --timestamping -P ./MM9/ \'ftp://hgdownload.cse.ucsc.edu/goldenPath/mm9/chromosomes/*' ";
				system "gunzip ./MM9/*.gz";
			} else {print "\n Well, then I can't really do my job. Goodbye. \n\n"; exit; }
		}
	} else {  
		#DOWNLOAD GENOME IF GENOME ISNT RECOGNIZED
		$filename = "./" . $genome . "/chr2.subst.fa";
		unless (-e $filename) { 
			print "\n Genome to isolate SNP loci not available, do you want to download the $genome Genome if it's available?   ";
			my $input = "";
			$input = <STDIN>;
			if ($input =~ m/Y/i){
				system "wget --timestamping -P ./$genome/ \'ftp://hgdownload.cse.ucsc.edu/goldenPath/$genome/chromosomes/*' ";
				system "gunzip ./$genome/*.gz";
			} else {print "\n Well, then I can't really do my job. Goodbye. \n\n"; exit; }
		}
	}
	
	
	
	
	#DELETE OLD SEQUENCES FILE AND CREATE NEW FILE
	if (-e './primer3/sequences.txt') {unlink('./primer3/sequences.txt')};
	open(SEQUENCES, ">./primer3/sequences.txt");
	
	
	
	#EXTRACT 8000bp REGIONS FLANKING SNPS 
	print "\n Extracting DNA sequences surrounding your SNPs.... \n";
	
	unless (defined($testsnp_file)) { $testsnp_file = "TestSNPs.txt";}
	open(TESTSNPS, $testsnp_file);
	@testsnps = <TESTSNPS>;
	close TESTSNPS;
	
	print "@testsnps\n\n\n";
	
	
	
	#GROUP SNPS BY CHROMOSOME IN ORDER TO LOOKUP ALL WHILE CHROMOSOME IS IN RAM - Chromosomes Hash references chromosome arrays where SNPs are held
	%chromosome_refs = ();
	
	foreach $snpline (@testsnps){
		chomp $snpline;
		my @snp = split(" ", $snpline);
	
		my $chromosome = $snp[0];
		my $location = $snp[1];
	
		push ( @{$chromosome_refs{$chromosome}}, $location);	
	}
	
	#Foreach Chromosome, lookup each SNP
	foreach $chromosome (sort keys %chromosome_refs){
		print "\n\nWill now lookup the sequences of these SNP locations in $chromosome   :\n ";
		print  " @{$chromosome_refs{$chromosome}} " . "\n";
	
	
	uc $genome;
	
	#FORMAT CHROMOSOME FILE INTO ONE LINE - SEQUENCE ONLY - TO USE BINARY SEEK / UNPACK ROUTINE
		my $chromosomefile = "./" . $genome . "/" . $chromosome . ".fa";
		print "\n\nOpening $chromosomefile ....\n\n";
		
		$newchromosomefile = $chromosomefile; 

		unless ($chromosomefile =~ /formatted/ ){

			$outfile = $chromosomefile . '.formatted';
	
			if (Extract->is_fasta($chromosomefile) =~ /false/ ) {print "Genome file is not FASTA - exiting\n\n"; exit; }
	
			open (IN, $chromosomefile);
			close IN;
	
			Extract->format_genome($chromosomefile , $outfile);
	
			$newchromosomefile = $outfile;
		}

	

	
		foreach $location (@{$chromosome_refs{$chromosome}}){
		
			$seq_length = 10000;
	 
			$start = $location - 5000;
			$end = $start + $seq_length;
	
			#DEFINE LENGTH AND LOCATION AS BEGINNING AND END
		
			print "Isolating $seq_length bp sequence from $start to $end of $chromosome in file $newchromosomefile which is very long.......\n\n";
	
		    $sequence = Extract->get_sequence($newchromosomefile, $start, $seq_length);
		   
			#$sequence =~ tr/GCATgcat/ /;
			$sequence =~ tr/SYRWKMsyrwkm/N/;
	
			my $snpline = $chromosome . "\t" . $location;
	
			print "Isolated genome substring $seq_length bases long.\n\n";
			print SEQUENCES $snpline . "\n" . $sequence . "\n";
			print $sequence;
		}
	}
	
	close SEQUENCES;
	
	
	
	
	
	###OPENS PRIMER 3 OPTIONS AND LOADS INTO ARRAY
	open (PRIMER3OPTIONS, "<./primer3/snp_primer3_options.txt");
	@primer3options = <PRIMER3OPTIONS>;
	close PRIMER3OPTIONS;
	
	
	#LOADS SEQUENCES EXTRACTED FROM GENOMIC DATA
	open (SEQUENCES , "<./primer3/sequences.txt");
	@sequencesfile = <SEQUENCES>;
	close SEQUENCES;
	
	print "\nSequences array: would normally print sequences here\n"; #print @sequencesfile;
	
	$sequenceslines = scalar @sequencesfile;
	
	
	## REMOVE OLD PRIMERS ENDFILE IF IT EXISTS / PREPARE PRIMER3_OUTPUT FOR ENTRIES
	system "rm ./primer3/primer3_output.txt";
	open(PRIMER3_OUTPUT, ">>./primer3/primer3_output.txt");
	print PRIMER3_OUTPUT "Output from Primer Designer program using primer3 and bowtie";
	
	
	
	#FOREACH SEQUENCE, CONSTRUCT PRIMER3 OPTIONS FILE, RUN PRIMER3 AND SAVE OUTPUT TO 'PRIMER3_OUTPUT'
	
	for ($i = 0; $i < $sequenceslines; $i +=2) {
		chomp $sequencesfile[$i] ;
		my $id = $sequencesfile[$i];
		print "\nFinding primer for $id \n";
		my $sequence_id = "SEQUENCE_ID=" . $id . "\n";
		my $sequence_template = "SEQUENCE_TEMPLATE=" . $sequencesfile[$i+1];
	
		print "writing primer3 input file for $id \n\n";
		
		my @primer3_initial_options = @primer3options;
	
		open (OUT, ">./primer3/primer3_input");
		print OUT $sequence_id;
		print OUT $sequence_template;
		print OUT @primer3_initial_options;
		close OUT;
	
	
		#RUN PRIMER3 INSIDE THIS LOOP
		&run_primer3($id);
	}
	
	
	close PRIMER3_OUTPUT;
	
	&extract_primer_sequences;
	
	
	#################################### EXIT STRATEGY #################################
	print "\n\nGoodbye, So Long, Farewell Auf Weidersein\n\n";
	exit;

}



sub run_primer3 {

	my $id = $_[0];
	print "Running Primer3 for $id ............\n\n";

PRIMER:
	my @primer = `./primer3/src/primer3_core ./primer3/primer3_input`;	
	
	#BEGIN TO LOOSEN THRESHOLDS IF PRIMER3 DOESNT FIND ANY PRIMERS - ADD PARAMETERS OTHER THAN Tm!
	if ($primer[($#primer - 1)] =~ /=0/) {
		print "\nReworking primer3 parameters\n\n";
		
		open (PRIMER3_INPUT , "<./primer3/primer3_input");
		@primer3_input = <PRIMER3_INPUT>;
		close PRIMER3_INPUT;
		
		#SHIFT Tms OUT ONE DEGREE EVERY ITERATION - make more responsive by only shifting troubling parameter
		$starting_primer_min_tm = substr ($primer3_input[9], -3, 2);
		$starting_primer_max_tm = substr ($primer3_input[11], -3, 2);
			$next_primer_min_tm = $starting_primer_min_tm -1;
			$next_primer_max_tm = $starting_primer_max_tm +1;
				print "Changing Primer_Min_Tm from $starting_primer_min_tm to $next_primer_min_tm \n";
				print "Changing Primer_Max_Tm from $starting_primer_max_tm to $next_primer_max_tm \n";			
					$primer3_input[9] = "PRIMER_MIN_TM=" . $next_primer_min_tm . "\n";
					$primer3_input[11] = "PRIMER_MAX_TM=" . $next_primer_max_tm . "\n";	
				
		$starting_primer_product_min_tm = substr ($primer3_input[12], -3, 2);
		$starting_primer_product_max_tm = substr ($primer3_input[13], -3, 2);
			$next_primer_product_min_tm = $starting_primer_product_min_tm -1;
			$next_primer_product_max_tm = $starting_primer_product_max_tm +1;
				print "Changing Primer_Product_Min_Tm from $starting_primer_product_min_tm to $next_primer_product_min_tm \n";
				print "Changing Primer_Product_Max_Tm from $starting_primer_product_max_tm to $next_primer_product_max_tm \n";
            	    $primer3_input[12] = "PRIMER_PRODUCT_MIN_TM=" . $next_primer_product_min_tm . "\n"; 
            	    $primer3_input[13] = "PRIMER_PRODUCT_MAX_TM=" . $next_primer_product_max_tm . "\n";


		#REWRITE PRIMER3_INPUT FILE WITH NEW Tm PARAMETERS
		open (PRIMER3_INPUT , ">./primer3/primer3_input");
		print PRIMER3_INPUT @primer3_input;
		close PRIMER3_INPUT;
		
		#RETRY PRIMER3 WITH NEW PARAMETERS
		goto PRIMER;
	}	
	
	
	#IF PRIMERS ARE FOUND, PRIMER3 OUTPUT IS APPENDED TO PRIMER3_OUTPUT
	print "PRIMER3 OUTPUT for $id: \n\n";
	print PRIMER3_OUTPUT "\n\nPRIMER3 OUTPUT for $id: \n\n";
	print @primer;
	print PRIMER3_OUTPUT @primer;

	return @primer;
	
}




sub extract_primer_sequences {	
	#EXTRACT PRIMER SEQUENCES FROM PRIMER3 OUTPUT FILE and SAVES TO FASTA FORMATTED FILE
	open (PRIMER3_OUTPUT , "./primer3/primer3_output.txt");
	our @primer3_output = <PRIMER3_OUTPUT>;
	close PRIMER3_OUTPUT;
		
	system "rm ./final_primers.fa";
	open (PRIMER3_SEQUENCES , ">>./final_primers.fa");


	
	foreach $line (@primer3_output){
		chomp $line;
	
		if ($line =~ /^PRIMER3 OUTPUT for /){
		our $id = $line;
		$id =~ s/PRIMER3 OUTPUT for //;
		}
	
		if ($line =~ /^PRIMER_LEFT_0_SEQUENCE/){
		$left = $line;
		$left =~ s/PRIMER_LEFT_0_SEQUENCE=//;
		#print "Left Primer Sequence:\n" . $left . "\n";
		}
		
		if ($line =~ /^PRIMER_RIGHT_0_SEQUENCE/){
		$right = $line;
		$right =~ s/PRIMER_RIGHT_0_SEQUENCE=//;
		#print "Right Primer Sequence:\n" . $right . "\n";
		}
		
		if ($line =~ /^PRIMER_LEFT_0_PROBLEMS/){
			$left_problems = $line;
			$left_problems =~ s/PRIMER_LEFT_0_PROBLEMS=//;
		}
		
		if ($line =~ /^PRIMER_RIGHT_0_PROBLEMS/){
			$right_problems = $line;
			$right_problems =~ s/PRIMER_RIGHT_0_PROBLEMS=//;
		}
		
		if (defined($left) && defined($right) && defined($id)){
		
			print ">Left_Primer_0_for: $id\n" . $left . "\n";
			print ">Right_Primer_0_for: $id\n" . $right . "\n";
			
			print PRIMER3_SEQUENCES ">Left_Primer_0_for_id: $id\n"  . $left . "\n";
			print PRIMER3_SEQUENCES ">Right_Primer_0_for_id: $id\n"  . $right . "\n";
		
			if (defined($left_problems) | defined($right_problems)){
				print "\nMade the best primers possible, but there might still be some issues:\n";
			}
			if (defined($left_problems)) {
				print "\nPotential problems with left primer: " . $left_problems;
			}
			
			if (defined($right_problems)){
				print "\nPotential problems with right primer: " . $right_problems;
			}
		
			$left_primer = $left;
			$right_primer = $right;
			
			undef($left);
			undef($right);
			undef($id);
			undef($left_problems);
			undef($right_problems);
		}	
	}
	
	close PRIMER3_SEQUENCES;		
}


