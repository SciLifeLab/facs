#!/usr/bin/perl -w

# Builds Bloom filters with a given K-mer length, false positive rate
# and reference in a FASTA format
#
# Copyright 2009 Henrik Stranneheim


=pod

=head1 NAME

BloomBuild - Create a filter for LoadFACS

=head1 SYNOPSIS

  BloomBuild k fp referencegenome outfile

=head2 ARGUMENTS

=over 4

=item B<k>

The word size for K-mers inserted into the filter.

=item B<fp>

The false positive probability that you accept

=item B<referencefile>

A file containing a list of filenames for reference genome data.
Each line contains a filename. Each file is a Fasta file.

=item B<outfile>

Filename for the created Bloom filter.

=back

=head1 DESCRIPTION

Build Bloom filters from reference genomes in Fasta files. These
can later be queried using LoadFACS.

=head1 SEE ALSO

LoadFACS

=cut



use strict;
use Bloom::Faster;
use Pod::Usage;
use Pod::Text;


# Inputs K-mer length, false positive and references to create bloom filters
if (@ARGV != 4) {
  my $verbosity = 0;
  if (@ARGV == 0) {
    $verbosity = 2;
  }

  pod2usage({-message => "Not the right number of arguments.",
	     -verbose => $verbosity
	    });
}

my ($targetlength, $falseposratebloom, $Ref, $outfile) = @ARGV;

my $bloomfilter;

my $id;	
my $seq;
my %seqs;
my $headercount;

my $refid;
my @refs;
my @accessions;
my $refcount;

my $keycount;
my $filtercounter=0;

#Reads list of reference files
#CheckReference($Ref);

#Refseq($refs[$refid], $refid);
Refseq($Ref);

for $id (keys %seqs) {	
#  print STDERR "id = $id\n";
  #Adds keys to the filter
  Targetkeys($seqs{$id}, $refid); 

  $id .=$targetlength;

}

print STDERR "Accessions:\n\t", join("\n\t", @accessions),"\n";
print STDERR "Targetlength:","\t", $targetlength,"\n";
print STDERR "False Positive Rate Bloom Filter:","\t", $falseposratebloom,"\n"; 

#Saves Bloom filter 
Savebloom("$outfile$targetlength");	
print STDERR "Created filter file: \t$outfile$targetlength\n";


#Resets %seqs & undefines bloom filter
#Blank($refid);

##########
#End of main program
##########

##########
##Subroutines
# my $id - Capture header in reference file	
# my $seq - Capture sequence in reference file
# my %seqs - Stores header and sequence in reference file  
# my $headercount - Counts the number of headers in reference file
##########
# my $refid - Capture filename in reference list
# my @refs - Stores filename in reference list
# my $refcount - Counts the number of filenames in reference file
##########
# my $keycount - Number of keys in filter
# my $filtercounter - Counts the number of filters created
##########


sub Bloomfilter {
    
	$keycount = $_[0];
#	$$_[1] = Bloom::Faster->new({n => $keycount, e => $falseposratebloom});

	return Bloom::Faster->new({n => $keycount, e => $falseposratebloom});
    
}


sub Savebloom {
  my $fname = shift @_;
  $bloomfilter->to_file($fname);
    
}


sub Blank {

    %seqs=();
    undef($$_[1]);
}


sub Targetkeys {
    
  my ($seq) = $_[0]; 
    
	my $lengthseq = length($seq)-($targetlength-1); #Stops sliding window at the end

	for (my $i=0;$i<$lengthseq;$i++) {	

	$bloomfilter->add (substr ($seq, $i, $targetlength) );
	
    }

    print STDERR "Added reference keys to filter\n";
    return;
    
}


sub Refseq {

   open(FASTA, "<$_[0]") or die "Can't open $_[0]:$!, \n";

   $headercount = 0;
   @accessions = ();
    while(<FASTA>) {

	if (/>/) {
	  $seq ="";
	  chomp($id = $'); #'
	  push @accessions, $id;
	  $headercount++;
	}
	
	else {

	    chomp($seq .= $_);
	}
	
    }

    $seqs{$id}= $seq;
    print STDERR "Finished Reading Reference Sequence","\n";
    close(FASTA);
	#Determines bitvector length & creates vector
   $bloomfilter = Bloomfilter(length($seqs{$id}), $$_[1] );  
    return;
    
}

sub CheckReference {
  my $fn = shift;

  open(BLOOM, "<$fn") or die "Can't open $fn:$!, \n";

  while(<BLOOM>) {

    chomp($refid = $_);  #'
    $refcount++;

    push @refs, $refid;
	}

  close(BLOOM);
  if ($refcount == 0) {
    print STDERR "Did not find any files in reference file '$fn'\n";
  }

  return;
}


