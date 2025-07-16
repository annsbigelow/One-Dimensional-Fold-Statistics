#!/usr/bin/perl
use Getopt::Std;
use Sys::Hostname;

getopts("bcde:hmn:q:rs:tvw");

# Print help message if -h option is specified
if($opt_h) {
    print "Usage: pov-movie.pl <switches> <snapshot-directory> <header-number> [<param_filename>]\n\n";
    print "Switches:\n";
    print "-b             (Set transparent background)\n";
    print "-c             (Add cylinders)\n";
    print "-d             (Don't duplicate existing files)\n";
    print "-e <num>       (Only render every <num> frame)\n";
    print "-h             (Print this information)\n";
    print "-m             (Switch off the mesh)\n";
    print "-n <num>       (Render up to <num> threads)\n";
    print "-q <quality>   (Quality of rendering, 1=good, 3=extreme)\n";
    print "-r             (Render remotely in parallel)\n";
    print "-s <frame>     (Render a single frame)\n";
    print "-t             (Render mesh as triangles, without normals)\n";
    print "-v             (Verbose output)\n";
    print "-w             (Disable making a movie)\n";
    exit 0;
}

die "Two or three arguments required" unless @ARGV==2 || @ARGV==3;

# Set variables used for remote processing
if($opt_r) {
    $rdir="esim/mesh/fiber";
    open A,"../config/rhosts" or die "Can't open remote hosts\n";
    @nlist=();@nthr=();
    while(<A>) {
        next if /^#/;
        @c=split;
        if($#c>=1) {
            push @nlist,$c[0];
            push @nthr,$c[1];
        }
    }
    close A;
    $nodes=$#nlist+1;
    $queue=$nodes==1?1:0;$h=0;
}

# Set miscellaneous variables, and those used to control which frames are
# rendered
$dr=$ARGV[0];
$msh=$opt_t?"mtr":"msh";
$verb=$opt_v?"":">/dev/null 2>/dev/null";
$every=$opt_e?$opt_e:1;
$a=defined $opt_s?$opt_s:0;
$opt_n=$opt_s if defined $opt_s;

# Read the first line of the POV-Ray header file to get the rendering
# dimensions. Assemble the POV-Ray flags.
$opt_q=1 unless defined $opt_q;
die "POV quality out of range\n" if $opt_q<0 || $opt_q>3;
open B,"pov_headers/h".($ARGV[1]).".pov"
    or die "Can't open POV-Ray header file\n";
$_=<B>;
$tr=$opt_b?"+UA":"";
m/^\/\/ W=(\d*) H=(\d*)/ or die "Can't read rendering size\n";
$pov_opts="+W$1 +H$2 ".
          ("","+R3 +A0.01 -J","+R6 +A0.001 -J","+R9 +A0.0001 -J")[$opt_q];

# Read, process, and store the rest of the POV-Ray header file
$gpn=0;
while(<B>) {
    if($opt_c) {s/^CYL://;} else {next if /^CYL:/;}
    if($opt_m) {next if /^MSH:/;} else {s/^MSH://;}
    $gp[$gpn++]=$_;
}
close B;$gpn--;

# Loop over the available frames
while(-e "$dr/pts.$a") {

    # Assemble output filename and check for skipping/termination conditions
    $fn=sprintf "fr_%04d.png",$a;
    last if defined $opt_n && $a>$opt_n;
    $a++,next if defined $opt_d && -e "$dr/$fn";

    # Assemble the POV file
    $tf="rtemp$h.pov";
    $par=(@ARGV==3?$ARGV[2]:"-");
    open A,$opt_r?"| bzip2 -9 -c >$dr/$tf.bz2":">$dr/$tf" or die "Can't open temporary POV file\n";
    foreach $i (0..$gpn) {
        $_=$gp[$i];

        if(/^#include "msh.pov"/) {
            print A `./unpack $msh $dr $a -`;
        } elsif(/^#include "cyl.pov"/) {
            print A `./unpack $cyl $dr $a -`;
        } else {
            print A;
        }
    }
    close A;

    # Render the frame, forking jobs to remote processors if the "-r" option is
    # given
    $pov_cmd="nice -n 19 povray $tf -D $tr +O$fn $pov_opts";
    if($opt_r) {

        # Send the POV-Ray file to a node for processing
        $hst=$nlist[$h];
        print "Frame $a to $hst\n";
        exec "rsync -q $dr/$tf.bz2 $hst:$rdir; ".
             "ssh $hst \"cd $rdir; bunzip2 -f $tf.bz2 ; $pov_cmd \" $verb ; ".
             "rsync -q $hst:$rdir/$fn $dr ; ssh $hst \"rm $rdir/$fn $rdir/$tf\" " if ($pid[$h]=fork)==0;

        # Wait for one of the forked jobs to finish
        if ($queue) {
            $piddone=wait;$h=0;
            $h++ while $piddone!=$pid[$h]&&$h<$nodes;
            die "PID return error!\n" if $h>=$nodes;
        } else {
            $h++;$queue=1 if $h>=$nodes-1;
        }
    } else {

        # Run POV-Ray locally
        print "Frame $a\n";
        system "cd $dr; $pov_cmd $verb";
    }
    $a+=$every;
}

if($opt_r) {wait foreach 0..($queue?$nodes-1:$h);}

# Additional code to automatically make a movie
unless ($opt_w) {
    ($mf=$dr)=~s/\.odr//;
    unlink "$mf.mov";
    system "ffmpeg -r 30 -y -i $dr/fr_%4d.png -preset veryslow -c:v libx265 -crf 17 -pix_fmt yuv420p -tag:v hvc1 -movflags faststart $mf.mov";
}

# Alternative movie-making command for older H.264 standard (more compatible, less efficient)
#ffmpeg -r 30 -y -i fr%4d.png -preset veryslow -c:v libx264 -crf 17 -pix_fmt yuv420p -movflags faststart plain_demo.mov
