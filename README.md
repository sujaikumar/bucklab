# bucklab

Follow these steps to get annotated chimeras from raw fastq sequencing files

## Trim fastq and convert to uniqfasta format

The format that we want is
```
>nId-nCount-nLen
ATTCTCTTATATATCGGGATA
```
where
- nID is a numeric ID (typically 1 for the most abundant unique sequence, but could be any number)
- nCount is the number of times that unique sequence is seen
- nLen is the length of the unique sequence
```
>1-240655-23
TAGCTTATCAGACTGATGTTGAC
>2-161403-22
AGGCAAGATGCTGGCATAGCTG
>3-158167-77
CTTCGTGATCGATGTGGTGACGTCGTGCTCTCCCGGGCCGGGTCCGAGCCGCGACGGGCGAGGGGCGGACGTTCGTG
```
Example command:
```
cutadapt --cores 8 --max-n=0 --minimum-length=18 -q 20 \
  -a TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC sample1.fq.gz 2>sample1.cutadapt.log \
  | fastq_to_uniqfa.pl >sample1.cutadapt.uniq.fa'
```

## Find chimeras

The goal is to identify reads where the first part hits a mature miRNA sequence (or some other small RNA / query RNA category, we expect these to be the first part of the chimera), and the second part hits the target genome. Because the second (i.e target) part can be quite short at times and can map to multiple locations in the genome, I use ShortStack's `--mmap u` mode to assign multi mapping reads.

Steps in the script:

1. Find the reads where the first part matches a miRNA and the second part is long enough to keep:
    1. Use bowtie2 to map reads to the smallRNA database (eg mature miRNA). The bowtie2 index with that prefix should already exist.
    1. Select only those reads where the first part matches the mature miRNA, and the remaining part is longer than a minimum threshold (default 18bp)
    1. Output the second part of these reads (with modified headers) to a temporary fasta file called btmpXXXXXXX (where bXXXXXXX is randomly generated)
1. Search the second part of the original reads against the target genome
    1. Before we do this, we have to first find WHICH part of the second part matches the target genome.
       This is because shortstack uses bowtie (not bowtie2) for mapping reads, and bowtie does not have a local-alignment mode.
       If the second part of the read is 80 bp long, but only bases 30-60 map to the target genome, then shortstack will not be able to find the match.
       Therefore the next step is to use bowtie2 (which has a local-alignment mode), to map btmpXXXXXXX against the target genome,
       pull out the part that matches the target genome, and store the output in a new temp fasta file called b2tmpXXXXXXX
    1. Before
    


We now need to search

1. With the first temp fasta file btmpXXXXXXX.fa as the input, use bowti2 the TARGET genome database

## Get annotation files 
