#!/usr/bin/perl -w

#Loads query in FASTA format, K-mer length and Bloom Filters. Interrogates queries against filters and classifies queries onto genomes. #The algorithm loops trough all queries for one filter at a time. 
#Copyright 2010 Henrik Stranneheim

=pod

=head1 NAME

facs - Filter reads of DNA

=head1 SYNOPSIS

  facs k bloomfilterlistfile queryfile outprefix

=head2 ARGUMENTS

=over 4

=item B<k>

The word size for K-mers used by the filter.

=item B<bloomfilterlistfile>

File containing a list of Bloom filters built using I<bloombuild>

=item B<queryfile>

Fasta file with reads to be filtered.

=item B<outprefix>

Prefix for output filenames.

=back

=head1 DESCRIPTION

facs interrogates queries against filters and classifies queries
onto genomes. The algorithm loops trough all queries for one filter
at a time.

Results are written to three files

=over 4

=item + Reads matching a reference

=item + Reads not matching a reference

=item + Reads with reverse complement matching a reference.

=back

=head1 SEE ALSO

bloombuild

=cut


use strict;
use Bloom::Faster;
use Pod::Usage;
use Pod::Text;

if (@ARGV != 4) {
  my $verbosity = 0;
  if (@ARGV == 0) {
    $verbosity = 2;
  }

  pod2usage({-message => "Not the right number of arguments.",
	     -verbose => $verbosity
	    });
}



# Inputs K-mer length, reference list and queries parameters
my ($targetlength, $reffilterlist, $queryfile, $outprefix) = @ARGV;
my $queryheadercount=0;

Querysseqs($queryfile); #reads query sequences and stores it in %querys

Reffilterlist($reffilterlist); #reads reference genome list and stores it in @refs

my $refid;
my $refcount;
my @refs;

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
my @matchesid;
my @matchscore;
my $norevcheck=0;
my $n_shorts = 0; 
my $trackmatch=0;
my $trackposition=0;


for (my $refid=0;$refid<scalar(@refs);$refid++) { 

	#Loads Bloom filter from reference array
	LoadBloom($refs[$refid]);
	push @matchesid, $refs[$refid];  #Stores reference id for classification stats later

	print STDERR "Filter:","\t", $refs[$refid],"\n";
	print STDERR "Targetlength:","\t", $targetlength,"\n"; 
    	
   	for $queryid (keys %querys) {  


		if ( length($querys{$queryid})>60 ) { #Minimum length cut-off
	
		#Divides querie into K-mers, checks filter, calculates match score and classifies sequences
		Sortquerykeys($querys{$queryid}, $targetlength,$queryid);
		
			if ($norevcheck eq 0) {
			#Translates and reverses queries
			$querys{$queryid} =~tr/ATGC/TACG/;
			$querys{$queryid} = reverse $querys{$queryid};
			Sortquerykeys($querys{$queryid}, $targetlength,$queryid);
			}
		$norevcheck=0; 

		}
		else {
		$allshort{$queryid}= $querys{$queryid}; #Saves all queries that did not surpass length cut-off
		delete $querys{$queryid}; #Deletes queries that did not surpass length cut-off from further querying
		$n_shorts++;		
		}

    	}	
    
    	print STDERR "Finished with Classification of Query Keys","\n";
		
	#Writes matches to file
	Write_matches($outprefix . "_" . $refs[$refid]);	
	Print();
	#Resets parameters
	Blank();
		if ( $refid>=(scalar(@refs)-1) ) {   # prints all Mismatches to file

			for my $nomatchid (keys %querys) { #Collects all unclassified sequences 

			$allshort{$nomatchid} = $querys{$nomatchid};
	    
			}

		Write_mismatches($outprefix . "_" . $refs[$refid]);

		}
}


#Prints matches and match score
for (my $matchid=0;$matchid<scalar(@matchesid);$matchid++) { #repeats until no more file id

print $matchesid[$matchid], "\t"; 
 $matchid++;
   if ($matchesid[$matchid]) {
	print $matchesid[$matchid], "\n"; 
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
# my @matchesid - Prints all reference headers and number of classified queries
# my $loopcount=0 - Loop for classifiying both forward and complementary reverse
# my $n_shorts = 0 - Counts all short sequences 
# my $trackmatch=0 - Tracks when to give a K-mer length score
# my $trackposition=0 - Tracks position within length of K-mer within last match 
##########


sub LoadBloom {
 
 my $filtername = shift;
 $filter = new Bloom::Faster($filtername);

}


sub Print {

 print STDERR "Number of sequences in original query file:","\t", $queryheadercount,"\n";
 print STDERR "Number of remaining queries:\t", scalar (keys %querys),"\n";
 print STDERR "Number of Mapped Targets:","\t", $bloommatchcount, "\n";
 print STDERR "Number of short queries: \t", $n_shorts, "\n";
 push @matchesid, $bloommatchcount;


}

sub Blank {

    $bloomnomatchcount =0;
    $bloommatchcount=0;
    $kmermatchcountqp=0;
    $kmermatchcount=0;
    $kmernomatchcount=0;
    @matches = ();
    @matchscore = ();     
}


sub Sortquerykeys {
    
    my $lengthseq = length($_[0]) -($_[1]-1); #Stops sliding window at the end
	
    for (my $i=0;$i<$lengthseq;$i++) {

	my $query = substr ($_[0], $i, $_[1]);

	#Quick pass of query
	CheckfilterQP($i, $query);

		if ($kmermatchcountqp) { 

	            if ($kmermatchcountqp eq 1) { #For queryies passing quickpass

			$kmermatchcountqp=0;
			$i=$lengthseq;
    			
			for (my $k=0;$k<$lengthseq;$k++) {

				my $query = substr ($_[0], $k, $_[1]);
	
				#Queries the Bloom filter
				Checkfilter($k, $query);

    			}
			
			if ( ( $kmermatchcount/length($_[0]) ) > 0.80) { #For matches
					
					push (@matches, $_[2]); #Saves match headers
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
	
	$_[0] = $_[0] + $targetlength-1; #Increments position in query 
	return;
    }
    
}

sub Checkfilter {

# $_[0] = Position in query
# $_[1] = my $query
     
    if ($trackposition >= $targetlength ) {

		$trackmatch = 0;
		$trackposition = 0;
    }

    if ($filter-> check( $_[1] ) ) {
	 
		if ($trackmatch eq 0) {
		
		$kmermatchcount = $kmermatchcount + $targetlength-1; #Adds to match score		
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
    if (m/^\s+$/) {		# Avoid blank lines
      next;
    }	
	if (/>/) {
	    
	    $queryseq ="";
	    chomp($queryid = $_); 
	    $queryheadercount++;
	}	
	
	else {
	    chomp($queryseq .= $_);

	}
	$querys{$queryid}=$queryseq;
	 
    }
    
    $queryheadercount = scalar(keys(%querys));
    close(QUERY);
    print STDERR "Finished Reading Query Sequences","\n";
    return;
}

sub Reffilterlist {
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
	$filename .= "_matching.fasta";
	
    open (GENOME, ">$filename") or die "Can't write to $filename: $!\n";
    
    my $assemblysseq;

	while (@matches) { 

		my $matchid = pop @matches;
		$assemblysseq .= $querys{$matchid};
	
		print GENOME $matchid,"\n"; 
	
		for (my $i=0;$i<(length($assemblysseq)/60);$i++) {
	    
		    print GENOME substr($assemblysseq,$i*60,60),"\n";
	    
		}
		$assemblysseq="";
		delete $querys{$matchid};
		
	}
     close (GENOME);
    return;
}

sub Write_mismatches {
    my $filename = shift;
    $filename .= "_mismatch.fasta";

    open (GENOME, ">$filename") or die "Can't write to $filename: $!\n";

    my $assemblysseq;
    
    foreach my $id (keys %allshort) {
	
	$assemblysseq .= $allshort{$id};
	
	print GENOME $id,"\n";
	
	for (my $i=0;$i<(length($assemblysseq)/60);$i++) {
	    
	    print GENOME substr($assemblysseq,$i*60,60),"\n";
	    
	}
	
	$assemblysseq="";
	
    }
     close (GENOME);
    return;
}

