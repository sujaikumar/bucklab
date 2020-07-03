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
eg:
```
>1-240655-23
TAGCTTATCAGACTGATGTTGAC
```
Example command:
```
cutadapt --cores 8 --max-n=0 --minimum-length=18 -q 20 \
  -a TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC sample1.fq.gz 2>sample1.cutadapt.log \
  | fastq_to_uniqfa.pl >sample1.cutadapt.uniq.fa'
```
