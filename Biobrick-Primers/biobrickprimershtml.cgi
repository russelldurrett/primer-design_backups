#!/usr/bin/perl 

#Thanks for taking a look at my code! - Russell Durrett ( red272 at nyu.edu ) 

#This script will take the input field 'DNA' from a form, calculate the 20bp sequence of forward and reverse primers and print them in HTML format.
#It can also calculate Tms, but since they are all 20bp it's not really necessary. 

#Feel free to cannablilize this code.


use CGI qw(:standard  escapeHTML);

my $DNA = param('DNA');


#Defines forward primer as first 20 bases

my $fprimer = substr ($DNA, 0, 20);
	
#Defines reverse primer as last 20 bases

&rprimer;
my $rprimer = substr ($DNArc, 1, 20);

sub rprimer {
	$_ = reverse $DNA;
	tr /AGCTagct/TCGAtcga/d;
	$DNArc = $_;
}


$fext = "CGCGCGAGAATTCGCGGCCGCTTCTAGA";
$rext = "GCGCGCACTGCAGCGGCCGCTACTAGT";


my $Tmfprimer = &Tm($fprimer);
my $Tmrprimer = &Tm($rprimer);

my $fprimer_wext = $fext . $fprimer;
my $rprimer_wext = $rext . $rprimer;


print "Content-type: text/html\r\n\r\n";
print  start_html("Biobricking Primers"), p("So your DNA was: ", tt(escapeHTML($DNA))), p("Your forward primer's sequence is: "), p(tt(escapeHTML($fprimer_wext))), p("Your reverse primer's sequence is: "), p(tt(escapeHTML($rprimer_wext))),  p("Primers are made with a 20bp overlap, thus all Tms should be above 50 degrees. If you really want to calculate the Tm, then download the source code from www.russelldurrett.com/biobrickprimers.txt and run it on your end."),  end_html();
exit;


##Perl script to calculate Tm of oligos - Not executed in this version. 

## Nearest Neighbor Tm subroutines were courtesy of a script written by David Gresham

sub Tm { 

my $seq = shift @_; 
 
$seq = uc($seq);

#print &calcTm(&calcNnEnthalpy($seq), &calcNnEntropy($seq));

}

#######subroutine to calculate Tm based on total enthalpy and entropy

sub calcTm {

    my ($totalEnthalpy, $totalEntropy) = @_;

    my $R = 1.9872; #gas constant cal/K-mol
    my $x = 4; #equals 4 for nonself-complementary duplex
    my $Ct = 0.0000000000006; #$Ct is strand concentration 
    
    #return ($totalEnthalpy * 1000) / ($totalEntropy + $R * log($Ct/$x)) - 273.15;

   
}



#######subroutine to calculate enthalpy by summing over nearest neighbors using the unified nearest neighbor values                  
sub calcNnEnthalpy {

#print "Calculating Enthalpy..\n";

    my $dna = shift(@_);

    my %enthalpies = (
        'AA' => -7.9,
        'TT' => -7.9,                                                                                                        
        'AT' => -7.2,                                                                                                        
        'TA' => -7.2,                                                                                                        
        'CA' => -8.5,                                                                                                        
        'TG' => -8.5,                                                                                                        
        'GT' => -8.4,                                                                                                        
        'AC' => -8.4,                                                                                                        
        'CT' => -7.8,                                                                                                        
        'AG' => -7.8,                                                                                                        
        'GA' => -8.2,                                                                                                        
        'TC' => -8.2,                                                                                                        
        'CG' => -10.6,                                                                                                      
        'GC' => -9.8,                                                                                                       
        'GG' => -8.0,                                                                                                        
        'CC' => -8.0                                                                                                          
        );

#select initialization value depending on whether terminal bases are A/T                                                            
#initialization and terminal penalty are not relevant for array hybridization 

    my $totalEnthalpy;

    if((substr($dna,0,1) eq 'A' || substr($dna,0,1) eq 'T') && (substr($dna,-1) eq 'A'|| substr($dna,-1) eq 'T'))
    {$totalEnthalpy = 0;} # 4.6;}                                                                                                   

    elsif((substr($dna,-1) eq 'A'|| substr($dna,-1) eq 'T') && (substr($dna,-1) ne 'A'|| substr($dna,-1) ne 'T'))
    {$totalEnthalpy = 0;} # 2.4;}                                                                                                   

    elsif((substr($dna,-1) ne 'A' || substr($dna,-1) ne 'T') && (substr($dna,-1) eq 'A' || substr($dna,-1) eq 'T'))
    {$totalEnthalpy = 0;} #2.4;}                                                                                                    

    else{$totalEnthalpy = 0;} #= 0.2};                                                                                              

#calculate total enthalpy based on each NN for oligo                                                                                

    for(my $i=0; $i< (length($dna) - 1); $i++) {

        $totalEnthalpy += $enthalpies{substr($dna,$i,2)};

    }

    return $totalEnthalpy;

}

#####subroutine to calculate entropy by summing over nearest neighbors using the nearest neighbor parameters
#initialization and terminal penalty are not relevant for array hybridization                         

sub calcNnEntropy {

#print "Calculating Entropy..\n";

    my $dna = shift(@_);

    my %entropies = (
        'AA' => -22.2,
        'TT' => -22.2,
        'AT' => -20.4,
        'TA' => -21.3,
        'CA' => -22.7,
        'TG' => -22.7,
        'GT' => -22.4,
        'AC' => -22.4,
        'CT' => -21.0,
        'AG' => -21.0,
        'GA' => -22.2,
        'TC' => -22.2,
        'CG' => -27.2,
        'GC' => -24.4,
        'GG' => -19.9,
        'CC' => -19.9
        );

#initialize total entropy depending on whether there are terminal A/T                                                                                    

    my $totalEntropy;

    if((substr($dna,0,1) eq 'A' || substr($dna,0,1) eq 'T') && (substr($dna,-1) eq 'A'|| substr($dna,-1) eq 'T'))
    {$totalEntropy = 0;} #8.1;}                                                                                                                          

    elsif((substr($dna,-1) eq 'A'|| substr($dna,-1) eq 'T') && (substr($dna,-1) ne 'A'|| substr($dna,-1) ne 'T'))
    {$totalEntropy = 0;} #1.2;}                                                                                                                          

    elsif((substr($dna,-1) ne 'A' || substr($dna,-1) ne 'T') && (substr($dna,-1) eq 'A' || substr($dna,-1) eq 'T'))
    {$totalEntropy = 0;} #1.2;}                                                                                                                          

    else {$totalEntropy = 0;} #-5.7;}                                                                                                                    

#add together entropies for each NN in string                                                                                                            

    for(my $i=0; $i < (length($dna) - 1); $i++) {

        $totalEntropy += $entropies{substr($dna,$i,2)};
    }

    return $totalEntropy;

}
#######subroutine to calculate deltaG at 37C by summing over nearest neighbors using the unified nearest neighbors values

sub calcDeltaG37 {

    my $dna = shift(@_);

    my %deltaG37 = (
        'AA' => -1.0,
        'TT' => -1.0,
        'AT' => -0.88,
        'TA' => -0.58,
        'CA' => -1.45,
        'TG' => -1.45,
        'GT' => -1.44,
        'AC' => -1.44,
        'CT' => -1.28,
        'AG' => -1.28,
        'GA' => -1.30,
        'TC' => -1.30,
        'CG' => -2.17,
        'GC' => -2.24,
        'GG' => -1.84,
        'CC' => -1.84
        );

    my $totalDeltaG37;

    if((substr($dna,0,1) eq 'A' || substr($dna,0,1) eq 'T') && (substr($dna,-1) eq 'A'|| substr($dna,-1) eq 'T'))
    {$totalDeltaG37 = 2.06;}

    elsif((substr($dna,-1) eq 'A'|| substr($dna,-1) eq 'T') && (substr($dna,-1) ne 'A'|| substr($dna,-1) ne 'T'))
    {$totalDeltaG37 = 2.01;}

    elsif((substr($dna,-1) ne 'A' || substr($dna,-1) ne 'T') && (substr($dna,-1) eq 'A' || substr($dna,-1) eq 'T'))
    {$totalDeltaG37 = 2.01;}

    else {$totalDeltaG37 = 1.96;}

    for(my $i=0; $i < (length($dna) - 1); $i++) {
  
        $totalDeltaG37 += $deltaG37{substr($dna,$i,2)};
    }

    return $totalDeltaG37;

}
