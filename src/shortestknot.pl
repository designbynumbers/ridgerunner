#!/usr/bin/perl

# This program scans a list of vect files and returns the shortest one.

use strict;
use warnings;
use File::Copy;
use File::Path;

my $myfile;
my $ropelength;
my $rbest = 10000000;
my $fbest = "";
my @files;

@files = <*.vect>;

print "Checking $#files files...";

my $rlcount = 0;

foreach $myfile (@files) {
    
    if ( $myfile =~ m/(\.vect)$/ && $myfile !~ m/(\.best\.vect)$/ 
	 && $myfile !~ m/(\.dVdt.vect)$/ && $myfile !~ m/(\.dlen.vect)$/ ) {

	$rlcount++;
    
	open(ROPELENGTH,"ropelength $myfile |");
	while (<ROPELENGTH>) {
	    if (/Ropelength:\s*(\S+)/) { 
		$ropelength = $1;
	    }
	}
	close ROPELENGTH;
	
	if ($ropelength < $rbest) {
	    
	    $rbest = $ropelength;
	    $fbest = $myfile;
	    
	}
    
    }

}

print "done.\n";
print "Computed ropelength for $rlcount files.\n";
print "Best ropelength value is $rbest for file $fbest.\n";

my $bestname = "";
my @fnparts;
my $npart;

@fnparts = split(/\./,$fbest);
pop(@fnparts);             # Remove the trailing .vect and the file number 
pop(@fnparts);
$bestname = join('.',@fnparts);

$bestname = $bestname.".best.vect";

copy($fbest,$bestname) or die("Could not copy $fbest onto $bestname.\n");
print("Copied best file $fbest to $bestname.\n");

