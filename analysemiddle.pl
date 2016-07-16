#!/bin/perl -w

use strict;
use Data::Dumper;

my $LEN = ($ARGV[0] or 100);

my %code;

while (<STDIN>) {
    chomp;
    if ($_ =~ /\"(.*)\"\s*,\s*(.*)\s*,\s*(.*)\s*,\s*(.*)\s*,\s*(.*)\s*,\s*(.*)\s*,\s*(.*)\s*,\s*(.*)\s*,\s*(.*)/) {
        $code{$1}->[$3] = [undef,$4,$5,$6,$7,$8,$9];
    }
}

my @stats;

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

for my $key (keys %code) {
    my @add = ($key);
    my $incomplete = 0;
    for my $e (1..6) {
        my $total = 0;
        my $sum = 0;
        for my $l ($e..$LEN) {
            if (not exists $code{$key}->[$l]) {
                $incomplete = 1;
            } else {
                my $weight = ($e == 1 ? ($l == 1 ? 1 : 0) : comb($e - 2, $l - 2)) * (31.0 ** $e) * ($LEN + 1 - $l);
                $sum += $code{$key}->[$l]->[$e] * $weight;
                $total += $weight;
            }
        }
        if (!$incomplete) {
            $add[$e] = $sum / $total;
        }
    }
    if (not $incomplete) {
        push @stats,\@add;
    }
}

@stats = sort { $a->[1] <=> $b->[1] or $a->[2] <=> $b->[2] or $a->[3] <=> $b->[3] or $a->[4] <=> $b->[4] or $a->[5] <=> $b->[5] or $a->[6] <=> $b->[6] } @stats;

print Dumper(\@stats);
