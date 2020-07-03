#!/usr/bin/env perl

=head1 NAME

mirfirst_chimeras.pl - Takes a uniq fa file from a clash chimera experiment, checks if the first part matches a miRNA
and returns a table with details about each chimera

=head1 SYNOPSIS

mirfirst_chimeras.pl \
  -i input_uniq.fa \
  -s small_rna_like_miRNA.fa \
  -t target_genome.fa \
  -n input.small.target.m18 \
  -m 18 \
  -b "--norc"
  
=head1 DESCRIPTION

=head1 AUTHORS

sujai.kumar@ed.ac.uk 2020-04-03

=cut

use strict;
use warnings;
use Getopt::Long qw(:config pass_through no_ignore_case);
use List::Util qw/ sum /;
use File::Temp qw/ tempfile tempdir /;

my $input_uniqfa_file;
my $srna_file;
my $target_genome_file;
my $name;
my $min_length = 18;
my $bowtie2first = "--norc";

GetOptions (
  "input_uniqfa_file=s" => \$input_uniqfa_file,
  "srna_file=s"         => \$srna_file,
  "target_genome_file=s"=> \$target_genome_file,
  "name=s"              => \$name,
  "min_length=s"        => \$min_length,
  "bowtie2first=s"      => \$bowtie2first,
);

die "-n name necessary" unless $name;

my $input_uniqfa_hash  = &fastafile_to_hash ($input_uniqfa_file);

# run bowtie2 against srna_file (bowtie2 index with that same name as prefix should exist)
#TODO: check bowtie2 index exists

open BOWTIE2, "bowtie2 --threads 8 -x $srna_file -f $input_uniqfa_file --no-unal --no-head -L 10 -i C,1 -N 1 --score-min C,30 --sensitive-local $bowtie2first |" or die $!;
my   ($fatmp_fh, $fatmp_filename) = tempfile("btmpXXXXXXX", DIR => ".");
my   %mirfirst;

while (<BOWTIE2>) {

  next if /^@/;
  my ($qname, $qstart, $qend, $qlen, $rname, $rstart, $rend, $rlen, $strand) = &extract_alignment_coordinates_from_sam ($_);
  
  # check enough bases left in read after srna_first:
  my $read_length = length($$input_uniqfa_hash{$qname}{seq});
  next if $read_length - $qend < $min_length;

  $mirfirst{$qname}{"srna_part_st" } = $qstart ;
  $mirfirst{$qname}{"srna_part_en" } = $qend ;
  $mirfirst{$qname}{"srna_ref_name"} = $rname ;
  $mirfirst{$qname}{"srna_ref_st"  } = $rstart ;
  $mirfirst{$qname}{"srna_ref_en"  } = $rend ;

  print $fatmp_fh ">${qname}:" . ($qend+1) . ":$read_length\n". substr($$input_uniqfa_hash{$qname}{seq}, $qend , $read_length - $qend) . "\n";
}

close $fatmp_fh;

close BOWTIE2;

# run second bowtie2 against target genome database to get the part of the read that needs to be shortstacked
# (bowtie2 index with that same name as prefix should exist)
#TODO: check bowtie2 index exists

open  BOWTIE2, "bowtie2 --threads 8 --reorder -x $target_genome_file -f $fatmp_filename --no-unal --no-head -L 10 -i C,1 -N 1 --score-min C,30 --sensitive-local |" or die $!;
my   ($fa2tmp_fh, $fa2tmp_filename) = tempfile("b2tmpXXXXXXX", DIR => ".");

my %max_as; # max alignment score
my ($as, $id, $orig_id, $orig_st, $orig_en, $new_id, $new_st, $new_en);
my ($qname, $qstart, $qend, $qlen, $rname, $rstart, $rend, $rlen, $strand);

while (<BOWTIE2>) {
  my $line = $_;
  next if $line =~ /^@/;
  #if ( /^(\S+)-\d+\t.*AS:i:(\d+)/ ) {
  if ( $line =~ /^(\S+).*AS:i:(\d+)/ ) {
    $id = $1;
    $as = $2
  } else {
    next
  };
  next if exists $max_as { $id } and $as < $max_as {$id};
  $max_as {$id} = $as;

  ($qname, $qstart, $qend, $qlen, $rname, $rstart, $rend, $rlen, $strand) = &extract_alignment_coordinates_from_sam ($line);
  
  # check enough bases in target_genome hit:
  next if $qend - $qstart < $min_length;

  # rename read to get coordinates of target_genome hit
  # eg  if 2345-34-67:24:67 matched from 3 to 40, then new read id would be 2345-34-67:26:63
  ($orig_id, $orig_st, $orig_en) = split ":",$id;
  $new_st = $orig_st + $qstart - 1;
  $new_en = $orig_st + $qend   - 1;
  $new_id = $orig_id . ":" . $new_st . ":" . $new_en;
  print $fa2tmp_fh ">$new_id\n". substr($$input_uniqfa_hash{$orig_id}{seq}, $new_st - 1 , $new_en - $new_st + 1) . "\n";
}

# remove low complexity using dustmasker, get first high complexity run of $min_length, rename id
my   ($fa3tmp_fh, $fa3tmp_filename) = tempfile("b3tmpXXXXXXX", DIR => ".");
open DUST, "dustmasker -in $fa2tmp_filename -outfmt fasta | seqtk seq | paste - - |" or die $!;
while (<DUST>) {
  chomp;
  my ($orig_id, $orig_st, $orig_en, $seq) = split /:|\t/, $_;
  my ($new_st, $new_len, $new_seq);
  my $is_dust = 1;
  while ( $seq =~ /([ATGCN]+)/g) {
    if (length($1) >= $min_length) {
      $is_dust = 0;
      $new_st  = length $`;
      $new_len = length $1;
      $new_seq = $1;
      last;
    }
  }
  next if $is_dust;
  print $fa3tmp_fh "$orig_id:" . ($orig_st + $new_st) . "-" . ($orig_st + $new_st + $new_len - 1) . "\n$new_seq\n"
}
close $fa3tmp_fh;

# convert to multifa and run shortstack
system "cat $fa3tmp_filename | uniqfa_to_multifa.pl > $name.multi.fa";
system "ShortStack --ranmax 'none' --pad 10 --dicermin 18 --dicermax 32 --nohp --bowtie_cores 30 --mmap u --bowtie_m all --sort_mem 4G " .
       "--readfile $name.multi.fa --genomefile $target_genome_file --outdir ShortStack_noranmax.$name.mall_mincov1";


exit;

### END MAIN ###

#############################################################################

sub read_fh {
  my $filename = shift @_;
  my $filehandle;
  if ($filename =~ /gz$/) {
    open $filehandle, "gunzip -dc $filename |" or die $!;
  }
  else {
    open $filehandle, "<$filename" or die $!;
  }
  return $filehandle;
}

#############################################################################

sub fastafile_to_hash {
  my $fastafile     = shift @_;
  my %sequences;
  my $fh = &read_fh($fastafile);
  my $seqid;
  while (<$fh>)
  {
    next if /^\s*$/ or /^#/;
    if (/^>(\S+)/) {
      $seqid = $1;
    }
    else {
      chomp($sequences{$seqid}{seq} .= $_  );
    }
  }
  return \%sequences;
}

#############################################################################

#get coordinates of qstart qend rstart rend using code from extract_alignment_coordinates_from_sam.pl  

sub extract_alignment_coordinates_from_sam {
  my $sam_line = shift @_;
  my @fields = split /\t/, $sam_line;
  next unless @fields >= 10;
  my $qname  = $fields[0];
  my $flag   = $fields[1];
  my $rname  = $fields[2];
  my $rstart = $fields[3];
  my $cigar  = $fields[5];
  my $qseq   = $fields[9];

  # Query
  my $qstart = 1;
  $_ = $cigar;
  s/^(\d+)[SH]/$qstart += $1/e;
  my $qlen = 0;
  $_ = $cigar;
  s/(\d+)[M=XI]/$qlen += $1/eg;
  my $qend = $qstart + $qlen - 1;

  # Reference
  my $rlen = 0;
  $_ = $cigar;
  s/(\d+)[M=XDN]/$rlen += $1/eg;
  my $rend = $rstart + $rlen - 1;

  # Strand
  my $strand = 1;
  if ($flag & 0x10) { $strand = -1 } else { $strand = 1 };
  return ($qname, $qstart, $qend, $qlen, $rname, $rstart, $rend, $rlen, $strand);
}

#############################################################################

sub revcomp {
  my $nuc = shift @_;
  $nuc =  uc($nuc);
  $nuc =~ tr /ATGCN/TACGN/;
  $nuc = scalar(reverse($nuc));
  return $nuc;
}
