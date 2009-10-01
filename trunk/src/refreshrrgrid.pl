#!/bin/perl

# This program refreshes the ridgerunner installations on the lab machines in 640.

use strict;
use warnings;
use File::Slurp;
use File::Path;
use FileHandle;
use IPC::Open2;

print "refreshrrgrid.pl is used to update the ridgerunner (etc) installations at UGA.\n";
print "(If you are not Jason Cantarella, you should not run this script.)\n";

# Get password from the user.

my $password;

print "Password:";
$password = <STDIN>;
chomp($password);

# Now go to each machine and do the installation

my @machines = qw( happy.math.uga.edu dopey.math.uga.edu sleepy.math.uga.edu sneezy.math.uga.edu );
my $hostname;

foreach $hostname (@machines) {

  my $output;
  my $input;
  
  print "-------------------------------\n";
  print "Updating machine $hostname...\n";
  print "-------------------------------\n";
  
  open2(*Reader, *Writer, "ssh $hostname \"cd plcurve; svn update; ./reconf; ./configure; make; sudo make install; cd ..; cd octrope; svn update; ./reconf; ./configure;  make; sudo make install; cd ..; cd tsnnls; svn update; ./configure; ./reconf; make; sudo make install; cd ..; cd ridgerunner; svn update; ./reconf; ./configure; make; sudo make install; \" ");
  
  print Writer "$password\n";

  while ( <Reader> ) {

    print $_;

    if ($_ =~ m/Password\:/) {

      print Writer "$password\n";

    }

  }

}


