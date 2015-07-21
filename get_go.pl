#!/usr/bin/env /proj/sot/ska/bin/perl

use warnings;

# GOES proton monitoring system for Chandra.

# Robert Cameron
# April 2001
# Brad Spitzbart
# March 2011
#  imported to CVS
#  currently using G13 as primary and G15 as secondary

use Time::JulianDay;
use PGPLOT;

$odir = "/data/mta4/proj/rac/ops/GOES";

#pgbegin (0,'/xs',1,3);
pgbegin (0,"$odir/pgplot.gif/vgif",1,2);
#pgsch(1.2);
pgsch(1.5);

$froot = 'ftp://ftp.swpc.noaa.gov/pub/lists/pchan/';
$fn = 'pchan_5m.txt';

@mon = qw(Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec);

#%id = ('GOES-15' => 'Gs', 'GOES-14' => 'G14', 'GOES-13' => 'Gp');
%id = ('GOES-15' => 'Gs', 'GOES-13' => 'Gp');

if ($p = `/usr/bin/lynx -source $froot/Gs_$fn`) { open PF, ">$odir/G15$fn"; print PF $p };
#if ($p = `/opt/local/bin/lynx -source $froot/G14$fn`) { open PF, ">$odir/G14$fn"; print PF $p };
if ($p = `/usr/bin/lynx -source $froot/Gp_$fn`) { open PF, ">$odir/G13$fn"; print PF $p };
if ($p = `/usr/bin/lynx -source ftp://ftp.swpc.noaa.gov/pub/lists/particle/Gp_part_5m.txt`) { open PF, ">$odir/G13_part_5m.txt"; print PF $p };

$t0 = time();
$tstamp = `date`;
chomp $tstamp;

gplot('GOES-15');
#gplot('GOES-14');
gplot('GOES-13');

sub gplot {
    $sat = $_[0];
    $r = $id{$sat};
    $t = $t0;
    @pt = ();
    @AoA = ();
    foreach (1..3) {
	($d,$m,$y) = (gmtime($t))[3,4,5];
	$y += 1900;
	$m++;
	$f = sprintf "%4.4d%2.2d%2.2d_",$y,$m,$d;
	$f = $froot.$f.$r."_".$fn;
	$p = `/usr/bin/lynx -source $f`;
	next unless ($p);
	@p = split m!$/!, $p;
	@pp = grep { /^\d/} @p;
	unshift @pt,@pp;
	$t -= 86400;
    }
    die scalar(gmtime)." No proton data found for $sat from $froot.\n" unless (@pt);
    foreach (@pt) {
	@f = split;
	foreach (0..$#f) { push @{$AoA[$_]},$f[$_] };
    }
#   @yr  = @{$AoA[ 0]};
#   @mn  = @{$AoA[ 1]};
#   @dy  = @{$AoA[ 2]};
#   @hm  = @{$AoA[ 3]};
    @mjd = @{$AoA[ 4]};
    @sod = @{$AoA[ 5]};
    @p1  = map { ($_ > 0)? log($_)/log(10) : -5 } @{$AoA[ 6]};
    @p2  = map { ($_ > 0)? log($_)/log(10) : -5 } @{$AoA[ 7]};
#   @p2  = map { ($_ > 0)? log($_)/log(10) : -5 } @{$AoA[ 8]};
#   @p3  = map { ($_ > 0)? log($_)/log(10) : -5 } @{$AoA[ 9]};
    @p5  = map { ($_ > 0)? log($_)/log(10) : -5 } @{$AoA[10]};
#   @p6  = map { ($_ > 0)? log($_)/log(10) : -5 } @{$AoA[11]};
#   @p7  = map { ($_ > 0)? log($_)/log(10) : -5 } @{$AoA[12]};
#   @p8  = map { ($_ > 0)? log($_)/log(10) : -5 } @{$AoA[13]};
#   @p9  = map { ($_ > 0)? log($_)/log(10) : -5 } @{$AoA[14]};
#   @p10 = map { ($_ > 0)? log($_)/log(10) : -5 } @{$AoA[15]};
#   @p11 = map { ($_ > 0)? log($_)/log(10) : -5 } @{$AoA[16]};

    $npt = $#pt + 1;
    $mjd0 = $AoA[4][0];
    $jd0 = $mjd0 + 2400001.5;
    @xlab = ();
    for (0..3) { 
	($xtit, $m0, $d0) = inverse_julian_day($jd0++);
	push @xlab,"$mon[$m0-1] $d0";
    }

    $xtit = "Coordinated Universal Time";
    $ytit = "p/cm2-s-sr-MeV";
    $ptit = "$sat Proton Flux";

    pgsci(1);
    pgpage;
    pgvstd;
    pgswin(0, 3, -3, 4);
    pgbox('BCTS', 0.5, 12, 'BCLNSTV2', 0, 0);
#    pgenv (0, 3, -3, 4, 0, 20);
    pglab ($xtit, $ytit, $ptit);

    $p4gm = log(300.0 / 3.3)/log(10);
    $p41gm = log(8.47 / 12.0)/log(10);

    @ephlx = (0, 3);
    @ephl4 = ($p4gm, $p4gm);
    @ephl41 = ($p41gm, $p41gm);

    @dd = @mjd;
    foreach (0..$#sod) { $dd[$_] += $sod[$_]/86400. - $mjd0 };

#    pgtext (1.99, 4.1, $tstamp);
    pgtext (2.1, 4.1, $tstamp);
    for (0..3) { pgptxt ($_, -3.3, 0, 0.5, $xlab[$_]) }; 

    pgsci(6);
    pgline ($npt, \@dd, \@p1);
    pgtext ( 0.0, 4.1, 'P1: 0.8-4');
    pgsci(3);
    pgline ($npt, \@dd, \@p2);
    pgtext ( 0.4, 4.1, 'P2: 4-9');
    pgtext (3.02, $p4gm+0.1, 'P4GM');
    pgtext (3.02, $p4gm-0.2, 'Limit');
    pgsci(5);
    pgline ($npt, \@dd, \@p5);
    pgtext (0.75, 4.1, 'P5: 40-80');
    pgtext (3.02, $p41gm+0.1, 'P41GM');
    pgtext (3.02, $p41gm-0.2, 'Limit');

    pgsci(2);
    pgline (2, \@ephlx, \@ephl4);
    pgline (2, \@ephlx, \@ephl41);
}

pgend;
