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

for my $key (keys %code) {
    my @add = ($key);
    my $incomplete = 0;
    for my $e (1..6) {
        my $total = 0;
        my $sum = 0;
        for my $l (1..$LEN) {
            if (not exists $code{$key}->[$l]) {
                $incomplete = 1;
            } else {
                $sum += $code{$key}->[$l]->[$e] * ($LEN + 1 - $l);
                $total += $LEN + 1 - $l;
            }
        }
        $add[$e] = $sum / $total;
    }
    if (not $incomplete) {
        push @stats,\@add;
    }
}

@stats = sort { $a->[1] <=> $b->[1] or $a->[2] <=> $b->[2] or $a->[3] <=> $b->[3] or $a->[4] <=> $b->[4] or $a->[5] <=> $b->[5] or $a->[6] <=> $b->[6] } @stats;

print Dumper(\@stats);
