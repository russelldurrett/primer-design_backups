#!usr/bin/perl

# This script will construct cloning primers for any nucleotide sequence entered on the command line
# Requires the primer3 directory (and customized primer3 options file) to be in the script directory. 

# Usage: perl primers.pl AGTCACAGTGACCAGGTAGAACA



$sequence = $ARGV[0];

NEXT:

print "OK, we're designing primers for this sequence:\n";
print $sequence;
print "\n\n";



#FORMAT SEQUENCE INTO SOMETHING READABLE BY PRIMER3 (ONLY N's)
$sequence =~ tr/SYRWKMsyrwkm/N/;


###OPENS PRIMER 3 OPTIONS AND LOADS INTO ARRAY
open (PRIMER3OPTIONS, "<./primer3/primer3_options.txt");
@primer3options = <PRIMER3OPTIONS>;
close PRIMER3OPTIONS;


## REMOVE OLD PRIMERS ENDFILE IF IT EXISTS / PREPARE PRIMER3_OUTPUT FOR ENTRIES
system "rm ./primer3/primer3_output.txt";
open(PRIMER3_OUTPUT, ">>./primer3/primer3_output.txt");
print PRIMER3_OUTPUT "Output from Primer Designer program using primer3 and bowtie";


#CONSTRUCT PRIMER3 OPTIONS FILE, RUN PRIMER3 AND SAVE OUTPUT TO 'PRIMER3_OUTPUT'
my $product_size = length $sequence;
print "Product size here is :   $product_size\n\n";

my $sequence_template = "SEQUENCE_TEMPLATE=" . $sequence . "\n";
my $target = "SEQUENCE_INCLUDED_REGION=" . "0, " . $product_size . "\n";

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







################ SUBROUTINES ################


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
		
		if (defined($left) && defined($right) && defined($id)){
		
		print ">Left_Primer_0_for: $id\n" . $left . "\n";
		print ">Right_Primer_0_for: $id\n" . $right . "\n";
		
		print PRIMER3_SEQUENCES ">Left_Primer_0_for_id: $id\n"  . $left . "\n";
		print PRIMER3_SEQUENCES ">Right_Primer_0_for_id: $id\n"  . $right . "\n";
		
		$left_primer = $left;
		$right_primer = $right;
		
		undef($left);
		undef($right);
		undef($id);
		
		$snps++;		
		}	
	}
	
	close PRIMER3_SEQUENCES;	
	print "\nMade primers for $snps gene!\n\n";
	
}



sub add_biobrick_extensions{

	print "Adding Biobrick Extensions to each primer - RFC 23 compatible\n\n";

	$fext = "CGCGCGAGAATTCGCGGCCGCTTCTAGA";
	$rext = "GCGCGCACTGCAGCGGCCGCTACTAGT";

	$final_left_primer = $fext . $left_primer;
	$final_right_primer = $rext . $right_primer;
}

