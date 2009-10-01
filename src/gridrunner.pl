#!/usr/bin/perl
#
# gridrunner. Given a command line argument listing a collection of files, generates an Xgrid batch file 
# which runs all of the files in the current directory with the given command line arguments (for ridgerunner).
# Unlike ridgerunner gridrunner requires the filename argument to come first.
#
# gridrunner <files> <args>
#
# The output file is always sent to stdout.
#

use strict;
use warnings;
use File::Slurp;
use Data::Dumper;
use File::Copy;
use Cwd;

my @rrargs;
my @files;

my $xgrid;

$xgrid = "xgrid -h happy.math.uga.edu -p \"X-grid\@640\" ";

#
# The first goal is to parse the command line arguments into files and rr args.
#

foreach (@ARGV) {

  if ($_ =~ m/vect/) {

    push(@files,$_);

  } else {

    push(@rrargs,$_);

  }

}

print STDERR "gridrunner 1.0, generating a batch file to run ".scalar @files." VECT files with rr arguments @rrargs \n";

#
# We now generate the header for the Xgrid batch file.
#

my $sec;
my $min;
my $hour;
my $mday;
my $mon;
my $year;
my $wday;
my $yday;
my $isdst;

($sec,$min,$hour,$mday,$mon,$year,$wday,
$yday,$isdst)=localtime(time);

my $cwd;
$cwd = getcwd;
# We need to delete the leading / from the $cwd string 
$cwd =~ s/^\///;

print "{\n";
print "jobSpecification = {\n";
print "  name = \" gridrunner batch job of ".scalar @files." files from $cwd generated at ";
printf "%4d-%02d-%02d %02d:%02d:%02d",$year+1900,$mon+1,$mday,$hour,$min,$sec;
print " \"; \n";

my $ridgerunner;
$ridgerunner = `which ridgerunner`;
chomp($ridgerunner);

#
# We now loop over files and assemble the header giving us access to the files for the run.
#	  
	
print "  inputFiles = { \n";

foreach (@files) {

#
# Now we need to trick xgrid into providing us with the input data for the 
# program in handy XML form. We do so by calling xgrid on a trivial run 
# with this input and then harvesting the specification file.
#
  my $jobnum_string; 

  $jobnum_string = `$xgrid -job submit $ridgerunner -s 10 $_`;
  $jobnum_string =~ /jobIdentifier = (\d+)/;

  my $jobnum;

  $jobnum = $1;

  my $jobspec  = `$xgrid -job specification -id $jobnum`;
  my $input_file;

  #print STDERR $jobspec;

  $jobspec =~ tr/\n/\#/;

  $jobspec =~ /inputFiles = ( .+? \}\;)/;
  $input_file = $1;
  $input_file =~ s/^.+?\"/\"/;
  $input_file =~ tr/\#/\n/;

  print "    $input_file\n\n";

}

print "};\n\n";

print "  taskSpecifications = { \n";

foreach (@files) { 

  my $shortname;
  my @narray;
  @narray = split('/',$_);
  $shortname = $narray[-1];
  $shortname =~ s/\W/\_/g;

  print "  $shortname = { \n";
  
  print "    command = \"$ridgerunner\"; \n";

  print "    arguments = ( \n";

  my $arg;
  foreach $arg (@rrargs) {
    print "      \"$arg\", \n";
  }
  print "      \"$cwd/$_\"\n";

  print "      ); \n";

  print "    };\n";

}

print "  };\n";

print "};\n";
print "}\n\n";

#
# Having done this, we try to generate a GridStuffer metajob file. This is going to be called gridstuffer.txt
#

open(GRIDSTUFFER,">gridstuffer.txt") or die("Couldn't open gridstuffer.txt.\n");

print GRIDSTUFFER "-dirs $cwd\n";

foreach (@files) {

  print GRIDSTUFFER "$ridgerunner $_ @rrargs\n";

}


