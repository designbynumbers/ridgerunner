#!/usr/bin/perl
#
#
# Animate a directory of vectfiles using povsnap. Should be run in the directory.
#

use warnings;
use File::Copy;
use Term::ProgressBar;

my $nframes;
$nframes = 20*24; # 20 seconds at 24 fps

my @vectfiles;
@vectfiles = <*.vect>;

my $skip;
$skip = int(scalar @vectfiles / $nframes);
if ($skip < 1) { $skip = 1; }

print "Animating directory of ".scalar @vectfiles." files, which we will cut down to ".(scalar @vectfiles / $skip)." frames to reach desired 20s time.\n";

my $file;
my $tubename;

my $kpname;
my $vectname;
my $pngname;

my $skipcount;
my $framecount;

$skipcount = 0;
$framecount = 1;

mkdir("../movie/");

my $progress;
$progress = Term::ProgressBar->new( {count => $nframes, ETA => linear,} );
$progress->max_update_rate(1);

foreach $file ( @vectfiles ) {

  if ($skipcount == $skip) { 

    $progress->update($framecount);

    `tube -r 0.5 $file 2>1 1>/dev/null`;
    
    $tubename = $file;
    $tubename =~ s/vect/tube\.off/g;
    
    `orient -a 2 $tubename`;
    `povsnap -v $tubename 2>1 1>/dev/null`;
    
    $pngname = $tubename;
    $pngname =~ s/off/png/g;
    
    move($pngname,"../movie/$pngname");
    unlink($tubename);

    $skipcount = 0;
    $framecount++;

  } else {

    $skipcount++;

  }

}
 
 
  
