#!usr/bin/perl

#package with subroutines to support Chromosomal Sequence Extraction without reading in the whole chromosome fileno

package Extract;

sub is_fasta {
	my $self = shift;
	my $file = shift;
	open (IN, $file);
	$first_character = getc IN;
	close IN;
	
	print "Running IS_FASTA? Routine - First character of file is $first_character";
	if ($first_character =~ />/ ) { return 'true' } else { return 'false' } 
}

sub format_genome {
	my $self = shift;
	my $genome_infile = shift;
	my $genome_outfile = shift;
	my $addition = ".formatted";
	
	my $formatted_genome_name = $genome_infile . $addition;
	
	unless (-e $formatted_genome_name){
		print "Reformatting $genome_infile to use with seek-unpack method. Renaming it $genome_outfile\n\n";
		
		open (IN, $genome_infile);
		@lines = <IN>;
		close IN;
		
		open (OUT, ">$genome_outfile");
		
		foreach $line (@lines){
		if ($line =~ /^>/){ #do nothing with header 
		} else {
			chomp $line;
			push @newlines, $line;
			}
		}
		$sequence = join ("", @newlines);
		
		print OUT $sequence;
		print "\n\nPRINTED ALL OF $genome_infile CONCATENATED TO ONE LINE NOW IN $genome_outfile \n\n\n";	
	}
}




sub get_sequence{


	#USAGE = get_sequence( INFILE , START POSITION , LENGTH );
	my $self = shift;
	my $infile = shift;
	my $start = shift;
	my $length = shift;
	
	open (IN, "<$infile");
	
	print "\n Getting $length characters from $infile starting at $start\n\n";
	
	$uplen = "A" . $length; #format length & input as ASCII Text (just 'A' and then the number of text characters)
	
	seek(IN, $start, 0); # seek to the position in the file you want to start at - beware newlines
	
	$output = unpack($uplen, <IN>); #read in bytes from binary
	
	return $output;

}









1;
