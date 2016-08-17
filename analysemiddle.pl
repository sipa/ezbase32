#!/bin/perl -w

use strict;
use Data::Dumper;
use List::Util qw( min max sum );

$| = 1;

my $LEN = ($ARGV[0] or 100);

sub comb {
    my ($k,$n) = @_;
    if ($k < 0) {
        return 0;
    }
    if ($n - $k < $k) {
        $k = $n - $k;
    }
    my $num = 1;
    my $den = 1;
    for my $i (1..$k) {
        $num *= ($n + 1 - $i);
        $den *= $i;
    }
    return ($num / $den);
}

my %code;
my @stats;

my @weights;

for my $e (1..12) {
    my $total = 0;
    for my $l ($e..$LEN) {
        my $weight = ($e == 1 ? ($l == 1 ? 1 : 0) : comb($e - 2, $l - 2)) * ($LEN + 1 - $l);
        $weights[$e]->[$l] = $weight / comb($e, $LEN);
    }
}

sub process {
    for my $key (keys %code) {
        my @add = ($key);
        my $incomplete = 0;
        for my $e (1..6) {
            my $total = 0;
            my $avg = 0;
            my $max = 0;
            for my $l ($e..$LEN) {
                if (not exists $code{$key}->[$l]) {
                    $incomplete = 1;
                } else {
                    if ($code{$key}->[$l]->[$e] > $max) {
                        $max = $code{$key}->[$l]->[$e];
                    }
                    $avg += $code{$key}->[$l]->[$e] * $weights[$e]->[$l];
                }
            }
            if (!$incomplete) {
                $add[$e] = [$avg, $max];
            }
        }
        if (not $incomplete) {
            push @stats,\@add;
        }
        delete $code{$key};
    }
}

my $key = "";
while (<STDIN>) {
    chomp;
    if ($_ =~ /\"(.*)\"\s*,\s*(.*)\s*,\s*(.*)\s*,\s*(.*)\s*,\s*(.*)\s*,\s*(.*)\s*,\s*(.*)\s*,\s*(.*)\s*,\s*(.*)/) {
        my $newkey = $1;
        my $len = $3;
        my $add = [undef, $4, $5, $6, $7, $8, $9];
        if ($key ne $newkey) {
            process();
            $key = $newkey;
        }
        $code{$newkey}->[$len] = $add;
    }
}
process();

sub bitcount {
    my ($n) = @_;
    my $ret = 0;
    while ($n > 0) {
        $ret += ($n & 1);
        $n = $n >> 1;
    }
    return $ret;
}

my $KEY = 0;

for my $stat (@stats) {
    for my $k (0..1) {
        for my $i (1..6) {
            print $stat->[$i]->[$k], " ";
        }
    }
    if ($stat->[0] =~ /0x(.*) 0x(.*) 0x(.*) 0x(.*) 0x(.*)/) {
        my @code = (hex $1,hex $2,hex $3,hex $4,hex $5);
        my @xcode = map { (($_&1)*$code[0]) ^ ((($_>>1)&1)*$code[1]) ^ ((($_>>2)&1)*$code[2]) ^ ((($_>>3)&1)*$code[3]) ^ ((($_>>4)&1)*$code[4]) } (1..31);
        my $mag = sum (map { log($_)/log(2) } @code);
        my $maxdev = max (map { abs(15 - bitcount($_)) } @xcode);
        my $sumdev = sum (map { abs(15 - bitcount($_)) } @xcode);
        print $mag, " ", $maxdev, " ", $sumdev, " ";
    } else {
        exit;
    }
    print "\"",$stat->[0], "\"\n";
}
