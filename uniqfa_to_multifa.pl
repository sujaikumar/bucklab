#!/usr/bin/env perl

open FA,"paste - - |";
while (<FA>) {
  if ( /^(>\S*\d+-(\d+)-\d+\S*)\t(\S+)$/ ) {
    for $i (1..$2) { print $1. "-$i\n$3\n" }
  } else {
    print "NOTHING"
  }
}
