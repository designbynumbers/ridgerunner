#!/bin/perl

use strict;
use warnings;
use File::Slurp;
use File::Path;
use Data::Dumper;
use File::Copy;
use Cwd;

#
#  prepare_tarball creates a complete ridgerunner tarball for distribution, including all support libraries.
#

sub pkgversion {
  
  my $pkgversion;
  my($path) = @_;

  open(CONFIGAC, "<", $path."configure.ac") 
	or die "cannot open < ${path}/configure.ac: $!";

  while (<CONFIGAC>) {

    if (m/AC_INIT\(.*,(.*),.*\)/) {

      $pkgversion = $1;
      last;

    }

  }

  return $pkgversion;

}

sub collectlibrary {

  my ($libname) = @_;
  my $libpath = "../".$libname."/";
  my $libver = pkgversion($libpath);
  my $makeoutput;
  my $libtar;
  my $libdir;

  $makeoutput = `cd $libpath; make dist;`;
  $makeoutput =~ m/.* -c >(.*.tar.gz)/;
  $libtar = $1;

  copy($libpath.$libtar,getcwd."/".$libtar) or die("Couldn't copy $libtar to ".getcwd."/".$libtar);

  print "$libtar -> ".getcwd."/${libtar}";
  `tar -xzf $libtar`;
  
  $libdir = $libtar;
  $libdir =~ s/\.tar\.gz//g;
  print " -> $libdir\n";

  unlink(getcwd."/${libtar}");

  return $libdir;
}

print "preparetarball 1.0\n\n";
print "Preparing ridgerunner tarball.\n";
print "Step 1. Refresh subpackages.\n";

my @libraries = qw( plCurve octrope vecttools plsurf povrayutils );
my $lib;
my $libdir;

foreach $lib (@libraries) {

  $libdir = collectlibrary($lib);
  die();

}




