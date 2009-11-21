#!/usr/bin/perl
#
#
# Discover the residual of a given vectfile, passing all arguments to 
# ridgerunner.
#

use strict;
use File::Temp;
use File::Copy;
use File::Slurp qw( slurp ) ;
use warnings;

my @files;
my @rrargs;

if (scalar @ARGV == 0) {

    print "Residual.pl computes the ridgerunner residual of a VECT file.\n";
    print "Usage: residual.pl <file.vect> <--lambda=(stiffness) > \n\n";
    print "Note: Any valid ridgerunner argument may be added to the residual\n";
    print "      command line (though most will have no visible effect).\n";
    die();

}

foreach (@ARGV) {

  if ($_ =~ m/vect/) {

    push(@files,$_);

  } 

  push(@rrargs,$_);

}

if ((scalar @files) == 0) {

    die("residual: Must provide at least one vect file in @ARGV.\n");

}

print STDERR "residual 1.0.\n\tComputing residual with rr command line\n";
print STDERR "\tridgerunner -c -s 1 --NoPNGOutput @rrargs \n";

my $tempdir;
my $rrdir;

$tempdir = File::Temp->newdir();

#print $tempdir;

$rrdir = $files[0];
$rrdir =~ s/vect/rr/; # replace extension

#print $rrdirname;

unless (-e $files[0]) {

    die("residual: File $files[0] does not exist.\n");

}

copy($files[0],$tempdir."/".$files[0]) or die("residual: Could not copy $files[0] to $tempdir.$files[0]");

chdir($tempdir) or die("residual: Could not chdir to $tempdir.\n");

my @rroutput;
@rroutput = `ridgerunner -c -s 1 --NoPNGOutput @rrargs`;

# print "ridgerunner output: @rroutput\n";

# We now open the residual.dat logfile to check the final residual.

my $resfile = slurp ($tempdir."/".$rrdir."/logfiles/residual.dat");

$resfile =~ /1 (.+)/g;
my $resnumber = $1;

print STDERR "Residual: ";
print $resnumber;












