#!usr/bin/perl


# THIS SCRIPT ALLOWS YOU TO QUERY GENBANK AND MAKE CLONING PRIMERS FOR ANY NUCLEOTIDE SEQUENCE
# AUTOMATICALLY ADDS BIOBRICK EXTENSIONS TO PRIMERS TO ALLOW CLONING, BUT OUTPUTS PRIMERS WITHOUT AS WELL

# Uses Bioperl to request GenBank files and store in queryable objects. (STREAM_OBJ)

# QUERY IS BASED ON ORGANISM AND GENE FIELDS, BUT CAN EASILYY BE MADE TO SUPPORT OTHERS

use lib './Modules/';

use Bio::Seq;
use Bio::SeqIO;
use Bio::DB::GenBank;
use Bio::DB::Query::GenBank;



$command_line_input = shift @ARGV;


#if ($command_line_input =~ /A|G|C|T/i && length($command_line_input) > 25){
# 	$sequence = $command_line_input;
# 	&cloning_primers;
# } 

if ($command_line_input =~ /.txt|.fa|.fasta/i){
	$input_fasta_file = $command_line_input;
	&cloning_primers_from_fafile;
} else { 
	&print_blank_input; 
	&genbank_primers;
}


sub print_blank_input{
	print "\n\nNo sequence or SNPs identified - running GenBank Query Protocol\n\n";
}




sub cloning_primers {
	
	
	NEXT:
	
	print "OK, we're designing primers for this sequence:\n";
	print $sequence;
	print "\n\n";
	
	
	
	#FORMAT SEQUENCE INTO SOMETHING READABLE BY PRIMER3 (ONLY N's)
	$sequence =~ tr/SYRWKMsyrwkm/N/;
	
	
	###OPENS PRIMER 3 OPTIONS AND LOADS INTO ARRAY
	open (PRIMER3OPTIONS, "<./primer3/cloning_primer3_options.txt");
	@primer3options = <PRIMER3OPTIONS>;
	close PRIMER3OPTIONS;
	
	
	## REMOVE OLD PRIMERS ENDFILE IF IT EXISTS / PREPARE PRIMER3_OUTPUT FOR ENTRIES
	system "rm ./primer3/primer3_output.txt";
	open(PRIMER3_OUTPUT, ">>./primer3/primer3_output.txt");
	print PRIMER3_OUTPUT "Output from Primer Designer program using primer3 and bowtie";
	
	
	#CONSTRUCT PRIMER3 OPTIONS FILE, RUN PRIMER3 AND SAVE OUTPUT TO 'PRIMER3_OUTPUT'
	my $product_size = length $sequence;
	print "Product size here is :   $product_size\n\n";
	
	my $end_base = $product_size -1;
	
	my $sequence_template = "SEQUENCE_TEMPLATE=" . $sequence . "\n";
	my $target = "SEQUENCE_INCLUDED_REGION=" . "0, " . $end_base . "\n";
	
	print "writing primer3 input file for $gene \n\n";
		
		open (OUT, ">./primer3/primer3_input");
		print OUT $sequence_template;
		print OUT $target;
		print OUT @primer3options;
		close OUT;
		
	
	#RUN PRIMER3 WITH THE NEW INPUT FILE
	&run_primer3($gene);
		
		
	&extract_primer_sequences;
	
	&add_biobrick_extensions;
	
	
	
	#PRINT FINAL PRIMERS TO FILE
	system "rm ./final_primers.fa";
	open (FINAL_OUTPUT , ">>./final_primers.fa");
	
	print FINAL_OUTPUT "> Forward Biobricking Primer for $gene :\n";
	print FINAL_OUTPUT $final_left_primer . "\n";
	print FINAL_OUTPUT "> Reverse Biobricking Primer for $gene : \n";
	print FINAL_OUTPUT $final_right_primer . "\n";
		
	print  "> Forward Biobricking Primer for $gene :\n";
	print  $final_left_primer . "\n";
	print  "> Reverse Biobricking Primer for $gene : \n";
	print  $final_right_primer . "\n";
	
	close FINAL_OUTPUT;
		
	
	
	
	
	####### Exit Strategy #######
	print "\n\nExiting............\n\n";
	exit;

}  #end Cloning Primers


sub cloning_primers_from_fafile {

	print "\nLoading fastafile $input_fasta_file sequences\n\n";
	$infile = Bio::SeqIO->new(-file => "$input_fasta_file" , -format => 'Fasta');

    while ( my $current_seq = $infile->next_seq() ) {
    
	
			NEXT:
			
			print "OK, we're designing primers for this sequence:\n";
			print $current_seq->id;
			print "\n\n";
			
			$gene = $current_seq->id;
			
			#FORMAT SEQUENCE INTO SOMETHING READABLE BY PRIMER3 (ONLY N's)
			$sequence = $current_seq->seq;
			$sequence =~ tr/SYRWKMsyrwkm/N/;
			
			
			###OPENS PRIMER 3 OPTIONS AND LOADS INTO ARRAY
			open (PRIMER3OPTIONS, "<./primer3/cloning_primer3_options.txt");
			@primer3options = <PRIMER3OPTIONS>;
			close PRIMER3OPTIONS;
			
			
			## REMOVE OLD PRIMERS ENDFILE IF IT EXISTS / PREPARE PRIMER3_OUTPUT FOR ENTRIES
			system "rm ./primer3/primer3_output.txt";
			open(PRIMER3_OUTPUT, ">>./primer3/primer3_output.txt");
			print PRIMER3_OUTPUT "Output from Primer Designer program using primer3 and bowtie";
			
			
			#CONSTRUCT PRIMER3 OPTIONS FILE, RUN PRIMER3 AND SAVE OUTPUT TO 'PRIMER3_OUTPUT'
			my $product_size = length $sequence;
			print "Product size here is :   $product_size\n\n";
			
			my $end_base = $product_size -1;
			
			my $sequence_template = "SEQUENCE_TEMPLATE=" . $sequence . "\n";
			my $target = "SEQUENCE_INCLUDED_REGION=" . "0, " . $end_base . "\n";
			
			print "writing primer3 input file for $gene \n\n";
				
				open (OUT, ">./primer3/primer3_input");
				print OUT $sequence_template;
				print OUT $target;
				print OUT @primer3options;
				close OUT;
				
			
			#RUN PRIMER3 WITH THE NEW INPUT FILE
			&run_primer3($gene);
				
				
			&extract_primer_sequences;
			
			&add_biobrick_extensions;
			
			
			
			#PRINT FINAL PRIMERS TO FILE
			#system "rm ./final_primers.fa";
			open (FINAL_OUTPUT , ">>./final_primers.fa");
			
			print FINAL_OUTPUT ">For $gene \t";
			print FINAL_OUTPUT $final_left_primer . "\n";
			print FINAL_OUTPUT ">Rev $gene \t";
			print FINAL_OUTPUT $final_right_primer . "\n";
				
			print  "> Forward Biobricking Primer for $gene :\n";
			print  $final_left_primer . "\n";
			print  "> Reverse Biobricking Primer for $gene : \n";
			print  $final_right_primer . "\n";
			
			close FINAL_OUTPUT;
				
	
	}
	
	
	####### Exit Strategy #######
	print "\n\nExiting............\n\n";
	exit;




}


sub genbank_primers {
	#Routine to identify organism / write in gene name query
	&choose_orgn_and_gene;
	
	
	#Submit search to GenBank and prompt for approval of good entry
	&submit_gene_search;
	
	
	NEXT:
	
	#PRINT SELECTED GENE DESCRIPTION, ACCESSION NUMBER AND NUCLEOTIDE
	
	print "OK, we're officially working with this guy:\n";
	print $final_obj->desc . "  Accession: " . $final_obj->primary_id;
	print "\n\nand its sequence:\n\n";
	print $final_obj->seq;
	print "\n\n";
	
	
	
	#FORMAT SEQUENCE INTO SOMETHING READABLE BY PRIMER3 (ONLY N's)
	$sequence = $final_obj->seq;
	$sequence =~ tr/SYRWKMsyrwkm/N/;
	
	
	###OPEN PRIMER 3 OPTIONS AND LOAD INTO ARRAY
	open (PRIMER3OPTIONS, "<./primer3/cloning_primer3_options.txt");
	@primer3options = <PRIMER3OPTIONS>;
	close PRIMER3OPTIONS;
	
	
	## REMOVE OLD PRIMERS ENDFILE IF IT EXISTS / PREPARE PRIMER3_OUTPUT FOR ENTRIES
	system "rm ./primer3/primer3_output.txt";
	open(PRIMER3_OUTPUT, ">>./primer3/primer3_output.txt");
	print PRIMER3_OUTPUT "Output from Primer Designer program using primer3 and bowtie";
	
	
	#CONSTRUCT PRIMER3 OPTIONS FILE, RUN PRIMER3 AND SAVE OUTPUT TO 'PRIMER3_INPUT'
	my $product_size = length $sequence;
	print "Product size here is :   $product_size\n\n";
	
	my $sequence_id = "SEQUENCE_ID=" . $gene . "\n";
	my $sequence_template = "SEQUENCE_TEMPLATE=" . $sequence . "\n";
	my $target = "SEQUENCE_INCLUDED_REGION=" . "0, " . $product_size . "\n";
	
	print "writing primer3 input file for $gene \n\n";
		
		open (OUT, ">./primer3/primer3_input");
		print OUT $sequence_id;
		print OUT $sequence_template;
		print OUT $target;
		print OUT @primer3options;
		close OUT;
		
	
	
	#RUN PRIMER3 WITH THE CUSTOMIZED INPUT FILE
	&run_primer3($gene);
		
		
	&extract_primer_sequences;
	
	&add_biobrick_extensions;
	
	
	
	#PRINT FINAL PRIMERS TO FILE
	system "rm ./final_primers.fa";
	open (FINAL_OUTPUT , ">>./final_primers.fa");
	
	print FINAL_OUTPUT ">For $gene \t";
	print FINAL_OUTPUT $final_left_primer . "\n";
	print FINAL_OUTPUT ">Rev $gene \t";
	print FINAL_OUTPUT $final_right_primer . "\n";
		
	print  "> Forward Biobricking Primer for $gene :\n";
	print  $final_left_primer . "\n";
	print  "> Reverse Biobricking Primer for $gene : \n";
	print  $final_right_primer . "\n";
	
	close FINAL_OUTPUT;

	
	
	####### Exit Strategy #######
	print "\n\Auf wiedersehen............\n\n";
	exit;
	
}	#end Genbank Primers




################ SUBROUTINES ################

sub choose_orgn_and_gene{
	print "\n\nThis script will find gene sequences and make Biobrick-compatible cloning primers for any nucleotide sequence in GenBank.\n\nWhich organism you would like to query:";

	$organism = <STDIN>;
	chomp $organism;
	
#	my $input = <STDIN>;
#	if ($input =~ /1/){
#		$organism = 'Escherichia coli';
#	} elsif ($input =~ /2/){
#		$organism = 'Deinococcus radiodurans';
#	} else {print "Couldn't understand that one...goodbye. \n\n"; exit;}


	print "\nWhat gene do you want to search for? ";
	$gene = <STDIN>;
	chomp $gene;
}




sub submit_gene_search{

	CYCLE: 

	$query = $gene . "[TITL] AND ". $organism . "[organism]";
	
	$query_obj = Bio::DB::Query::GenBank->new(
		-query => $query, 
		-db => 'nucleotide' );
	
	$gb_obj = Bio::DB::GenBank->new;
	
	$stream_obj = $gb_obj->get_Stream_by_query($query_obj);

	$counter = 0;
	
		#LIST ALL FILES
	while ($seq_obj = $stream_obj->next_seq) {
		$counter++;
		
		if ($counter == 1){
			print "\n\nHere's a list of all the compatible GenBank files to choose from:\n";
		}
		print "\n" . $counter . " " . $seq_obj->desc;
	}
	
	$stream_obj = $gb_obj->get_Stream_by_query($query_obj);
	$counter = 0;
		
		#AND THEN CYCLE THROUGH THEM TO PICK ONE
	while ($seq_obj = $stream_obj->next_seq) {
		$counter++;
		
		if ($counter == 1){
			print "\n\n\nNow we'll cycle through them - when you see the one you want enter 'good'\n\n";
		}
			
		print "\n How does this one look to you? Enter 'good' or press enter\n\n";
		
		print "#" . $counter . ": " . $seq_obj->desc . $seq_obj->primary_id . "\n";
		
	
		my $input = <STDIN>;
	
		if ($input =~ /good/i) {
			$final_obj = $seq_obj;
			print "Saved object as final target....\n\n";
			goto NEXT;
		}
	print "\n\n";
	}
	
	if ($counter == 0){
		print "\n\nNO MATCHES FOUND - SORRY.\n\n";
		exit;
	}
	
	
	print "No More Matches, do you want to cycle through them again?\n\n";
	my $cycle_again = <STDIN>;
	if ($cycle_again =~ /Y/i){
		goto CYCLE;
	} else {print "OK, goodbye then. \n\n"; exit;}
	
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
		
	#system "rm ./final_primers.fa";
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
			
			#print PRIMER3_SEQUENCES ">Left_Primer_0_for_id: $id\n"  . $left . "\n";
			#print PRIMER3_SEQUENCES ">Right_Primer_0_for_id: $id\n"  . $right . "\n";
		
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



sub add_biobrick_extensions{

	print "\n\nAdding Biobrick Extensions to each primer - these are now RFC 23 compatible:\n\n";

	$fext = "CGATCGAGAATTCGCGGCCGCTTCTAGA";
	$rext = "GCTATGCACTGCAGCGGCCGCTACTAGT";

	$final_left_primer = $fext . $left_primer;
	$final_right_primer = $rext . $right_primer;
}

