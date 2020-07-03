#!/usr/bin/env perl

open FQ,"awk 'NR%4==2'|";
while (<FQ>) {
  chomp;
  $s{$_}++;
}

foreach (sort {$s{$b} <=> $s{$a}} keys %s) {
  print ">".++$i."-".$s{$_}."-".length($_)."\n$_\n"
}
