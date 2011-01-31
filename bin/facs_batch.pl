#!/usr/bin/perl -w

# Loads query in FASTA format, K-mer length and Bloom Filters. Interrogates
# queries against filters and classifies queries onto genomes. 
#The algorithm loops trough all queries for one filter at a time. 
# Copyright 2011 Henrik Stranneheim

=pod

=head1 NAME

facs - Filter reads of DNA

=head2 SYNOPSIS

facs queryfile -b bloomfilterlist

=head3 DESCRIPTION

Build Bloom filters from reference genomes in fasta files using bloombuild.pl. These
can later be queried using FACS.

This version of facs.pl is based on facs_2.pl, using the more accurate scoring
system.

facs interrogates queries against filters and classifies queries
onto genomes. The algorithm loops trough all queries for one filter
at a time.

Results are written to two files

1. Reads matching a reference

2. Reads not matching any reference

=head3 COMMANDS AND OPTIONS

-b/--bloomfilter A file containing a list of filenames for already created bloom filters. Each line contains the file name of the specific bloom filter. (defaults to STDIN)

-op/--outfileprefix Output prefix for the output files (defaults to "")

-osma/--outfilesuffix Output suffix for the matches output files

-osmi/--outfilesuffix Output suffix for the mismatches output files

-k/--kmer Word length for K-mers to be queried against the bloom filter (defaults to '21')

-f/--falseposprob The false positive probaility that you accept (defaults to '0.0005')

-dbpr/--databsepath The unix path to the bloom filters (defaults to STDIN)

-o/--outdirectory The unix path to the output directory (defaults to STDIN)

-mc/--matchcutoff The percent identity to classify a match (defaults to 80)

-lc/--lengthcutoff The minimum required read length (defaults to 60)

-qp/--quickpasses Number of quick passes that must match before vetting read (defaults to 1)

-bs/--batchsize Number of reads that are analysed in batch (defaults to 5000)

Note: You have to use the same false positive rate and K-mer length as was used when the bloom filters were created.

=head4 I/O

Input format (FASTA/Pearson)

Output format (FASTA/Pearson)

=head5 SEE ALSO

bloombuild

=cut

use strict;
use Bloom::Faster;
use Pod::Usage;
use Pod::Text;
use Getopt::Long;

use vars qw($USAGE);

BEGIN {
    $USAGE =
qq{facs.pl < queryfile.fa -b bloomfilter.obj > file_suffix.fasta
-k/--kmer Word length for K-mers to be queried against the bloom filter (defaults to '21')
-f/--falseposprob The false positive probaility that you accept (defaults to '0.0005')
-b/--bloomfilter A file containing a list of filenames for already created bloom filters. Each line contains the file name of the specific bloom filter. (defaults to STDIN)
-op/--outfileprefix Output prefix for the output files (defaults to "")
-osma/--outfilesuffix Output suffix for the matching output files (defaults to "_matching.fasta")
-osmi/--outfilesuffix Output suffix for the mismatching output files (defaults to "_mismatching.fasta")
-dbpr/--databsepath The unixpath to the bloom filters (defaults to STDIN)
-o/--outputdirectory The unixpath to the output directory (defaults to STDIN)
-mc/--matchcutoff The percent identity to classify a match (defaults to "0.8"~80%)
-lc/--lengthcutoff The minimum required read length (defaults to 60)
-qp/--quickpasses Number of quick passes that must match before vetting read (defaults to 1)
-bs/--batchsize Number of reads that are analysed in batch (defaults to 5000)
};

}
my ($rformat, $oformat, $outfileprefix, $matchoutfilesuffix, $mismatchoutfilesuffix, $kmerlength, 
$falseposratebloom, $dbpr, $od, $matchcutoff, $lengthcutoff, $nrqp, $bs, 
$help) = ('fasta','fasta', "", "_matching.fasta", "_mismatching.fasta", 21, 0.0005,".", ".", 0.8, 60, 1, 5000);
my ($queryfile, $bloomfilterlist);

GetOptions('k|kmerlength:i' => \$kmerlength,
'f|falseposratebloom:s' => \$falseposratebloom,
'b|bloomfilter:s' => \$bloomfilterlist,
'op|outfileprefix:s' => \$outfileprefix,
'osma|matchoutfilesuffix:s' => \$matchoutfilesuffix,
'osmi|mismatchoutfilesuffix:s' => \$mismatchoutfilesuffix,
'dbpr|databaseref:s'  => \$dbpr,
'o|outdirectory:s'  => \$od,
'mc|matchcutoff:s'  => \$matchcutoff,
'lc|lengthcutoff:i'  => \$lengthcutoff,
'qp|quickpasses:i'  => \$nrqp,
'bs|batchsize:i'  => \$bs,
'h|help' => \$help,
);

die $USAGE if( $help );

if ($queryfile) {

	}
else {

($queryfile) = @ARGV;

if (@ARGV != 1) {
  my $verbosity = 0;
  if (@ARGV == 0) {
    $verbosity = 2;
  }
  print"\n";
  pod2usage({-message => "Must supply a query file and list of bloom filter(s).\n",
	     -verbose => $verbosity
	    });
	}
}

my $refid;
my $refcount;
my @refs;
my %refs;

my $queryid;
my %querys;
my $queryseq;

my $filter;
my $kmermatchcountqp=0;
my $kmermatchcount=0;
my $kmernomatchcount=0;
my $bloomnomatchcount=0;
my $bloommatchcount=0;
my %allshort;
my @matches;
my $norevcheck=0;
my $n_shorts = 0; 
my $trackmatch=0;
my $trackposition=0;

my $queryheadercount=0;

Bloomfilterlist($bloomfilterlist); #reads reference genome list and stores it in @refs

Querysseqs($queryfile); #reads query sequences and stores it in %querys

sub Classify {

for (my $refid=0;$refid<scalar(@refs);$refid++) { 

	#Loads Bloom filter from reference array
	LoadBloom($refs[$refid]);
    	
   	for $queryid (keys %querys) {  


		if ( length($querys{$queryid})>$lengthcutoff ) { #Minimum length cut-off
	
		#Divides querie into K-mers, checks filter, calculates match score and classifies sequences
		Sortquerykeys($querys{$queryid}, $kmerlength,$queryid, $refs[$refid]);
		
			if ($norevcheck eq 0) {
			#Translates and reverses queries
			$querys{$queryid} =~tr/ATGC/TACG/;
			$querys{$queryid} = reverse $querys{$queryid};
			Sortquerykeys($querys{$queryid}, $kmerlength,$queryid, $refs[$refid]);
			}
		$norevcheck=0; 

		}
		else {
		$allshort{$queryid}= \$querys{$queryid}; #Saves all queries that did not surpass length cut-off
		$n_shorts++;		
		}

    }	
    
		
	#Writes matches to file
	Write_matches($outfileprefix . "_" . $refs[$refid]);	
	#Resets parameters
	Blank();
		if ( $refid>=(scalar(@refs)-1) ) {   

			for my $nomatchid (keys %querys) { #Collects all unclassified sequences 

			$allshort{$nomatchid} = \$querys{$nomatchid};
	    
			}

		Write_mismatches($outfileprefix . "_");		# prints all Mismatches to file

		}
	}
	return;
}

Print();

for $refid (keys %refs) {		#Prints matches and match score

print $refid, "\t"; 

   if ($refs{$refid}) {
   
	print $refs{$refid},"\n";
   }
   else {
	print "\n"; 
   }

}

undef($filter);

##########
#End of main program
##########

##########
###Subroutines:
# my $refid - Capture filename in reference list
# my @refs - Stores filename in reference list
# my $refcount - Counts the number of filenames in reference file
##########
# my $queryid - Capture each header in queryfile
# my %querys - Stores header and corresponding sequences in queryfile
# my $queryseq - Capture each sequences in queryfile
# my $queryheadercount=0 - Counts number of queries
##########
# my $filter - Bloom filter object
# my $kmermatchcountqp=0 - Number of matching K-mers per query in Quick pass
# my $kmermatchcount=0 - Number of matching nucleotides (K-mers*K-mer length) per queries, match score
# my $kmernomatchcount=0 - Number of mismatching K-mers per queries
# my $bloomnomatchcount=0 - Number of classifed queries
# my $bloommatchcount=0 - Number of John Does
# my %allshort - Stores all queries below length cut-off and at the end collects all John Does for printing
# my @matches - Stores query id for classification 
# my %refs - Prints all reference headers and number of classified queries
# my $loopcount=0 - Loop for classifiying both forward and complementary reverse
# my $n_shorts = 0 - Counts all short sequences 
# my $trackmatch=0 - Tracks when to give a K-mer length score
# my $trackposition=0 - Tracks position within length of K-mer within last match 
##########


sub LoadBloom {
 
 my $filtername = shift;
 $filter = new Bloom::Faster("$dbpr/$filtername");

}


sub Print {
 print STDERR "Finished with Classification of Query Keys","\n";
 print STDERR "Number of sequences in original query file:","\t", $queryheadercount,"\n";
 print STDERR "Number of short queries: \t", $n_shorts, "\n";


}

sub Blank {

    $bloomnomatchcount =0;
    $bloommatchcount=0;
    $kmermatchcountqp=0;
    $kmermatchcount=0;
    $kmernomatchcount=0;
    @matches = ();
    #@matchscore = ();     
}


sub Sortquerykeys {

#$_[0] = $queryseq
#$_[1] = $kmerlength
#$_[2] = $queryid
#$_[3] = $refid
    
    my $lengthseq = length($_[0]) -($_[1]-1); #Stops sliding window at the end
	
    for (my $i=0;$i<$lengthseq;$i++) {

	my $query = substr ($_[0], $i, $_[1]);

	#Quick pass of query
	CheckfilterQP($i, $query);

		if ($kmermatchcountqp) { 

	            if ($kmermatchcountqp eq $nrqp) { #For queryies passing quickpass

			$kmermatchcountqp=0;
			$i=$lengthseq;
    			
			for (my $k=0;$k<$lengthseq;$k++) {

				my $query = substr ($_[0], $k, $_[1]);
	
				#Queries the Bloom filter
				Checkfilter($k, $query);

    			}
			
			if ( ( $kmermatchcount/length($_[0]) ) > $matchcutoff) { #For matches
					
					push @matches, \$_[2]; #Saves match headers
					$refs{$_[3]}++;
					$norevcheck=1;
				    $bloommatchcount++;
	   				$kmermatchcount=0;
	   				$kmernomatchcount=0;
					$trackmatch = 0;
					$trackposition = 0;
				 	return;
	  		}
    		$kmermatchcount=0;
   			$kmernomatchcount=0;
			$trackmatch = 0;
			$trackposition = 0;
  			return;
		}
	}
   }

}

sub CheckfilterQP {

# $_[0] = Position in query
# $_[1] = my $query
     
    if ($filter-> check( $_[1] ) ) {
	 
	 ++$kmermatchcountqp; #Adds to match score 
	return;
    }
    
    else {
	
	$_[0] = $_[0] + $kmerlength-1; #Increments position in query 
	return;
    }
    
}

sub Checkfilter {

# $_[0] = Position in query
# $_[1] = my $query
     
    if ($trackposition >= $kmerlength ) {

		$trackmatch = 0;
		$trackposition = 0;
    }

    if ($filter-> check( $_[1] ) ) {
	 
		if ($trackmatch eq 0) {
		
		$kmermatchcount = $kmermatchcount + $kmerlength-1; #Adds to match score		
		$trackmatch = 1;
		$trackposition++;
		}
		else {
		$kmermatchcount++;
		$trackposition++;
		}
	return;
    }
    
    else {
	
	$trackposition++;

	return;
    }
    
}


sub Querysseqs {

    open(QUERY, "<$_[0]") or die "Can't open $_[0]:$!, \n";    

    while (<QUERY>) {
        chomp $_;

        if (m/^\s+$/) {		# Avoid blank lines
            next;
        }	

	    if (/>(\S+)/) {
	    
	 		if (scalar(keys(%querys)) eq $bs) {   #Batch size
	    	
	    		print $queryheadercount=$queryheadercount + scalar(keys(%querys)), "\n";
	    		Classify();
	    		%querys = ();

	    	}
	        $queryseq ="";
	        $queryid = $1; 	        
	        
	    }	
	    else {
	        $queryseq .= $_;
	    }

	    $querys{$queryid}=$queryseq;
    }
    
    $queryheadercount=$queryheadercount + scalar(keys(%querys));
    Classify(); #Catch the remainaing reads
    close(QUERY);
    print STDERR "Finished Reading Query Sequences","\n";
    return;
}

sub Bloomfilterlist {
 my $filename = shift;

 open(REF, "<$filename") or die "Can't open $filename:$!, \n";   

    $refcount = 0;
    while(<REF>) {

	
	if (/\S+/) {

      chomp;
      $refcount++;
      push @refs, $_;
    }

    }
    
    close(REF);
    print STDERR "Finished Reading Ref Filter List","\n"; 
    return;
}


sub Write_matches {
    
	my $filename = shift;
	$filename .= $matchoutfilesuffix;
	
    open (GENOME, ">>$od/$filename") or die "Can't write to $od/$filename: $!\n";

	while (@matches) { 

		my $matchid = pop @matches;
		
		print GENOME '>', $$matchid,"\n"; 

	    for (my $i=0;$i<(length($querys{$$matchid})/60);$i++) {
	    
	    	print GENOME substr($querys{$$matchid},$i*60,60),"\n";
		}
		
		delete $querys{$$matchid};
		
	}

     close (GENOME);
    return;
}

sub Write_mismatches {
    my $filename = shift;
    $filename .= $mismatchoutfilesuffix;

    open (GENOME, ">>$od/$filename") or die "Can't write to $od/$filename: $!\n";
   
    foreach my $id (keys %allshort) {
	
	print GENOME '>', $id,"\n";
	
	for (my $i=0;$i<(length(${$allshort{$id}})/60);$i++) {
	    
	    print GENOME substr(${$allshort{$id}},$i*60,60),"\n";
	    
	}
	
	delete $allshort{$id};
    }
     close (GENOME);
    return;
}

