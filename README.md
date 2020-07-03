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

Each sample uniqfa file will typically be in its own SAMPLE directory. I'm using this directory structure:
```
project
|__sample1
|  |_sample1.cutadapt.uniq.fa
|
|__sample2
|  |_sample2.cutadapt.uniq.fa
|
|__sample3
|  |_sample3.cutadapt.uniq.fa
|
|__sample4
   |_sample4.cutadapt.uniq.fa
```
## Find chimeras

The goal is to identify reads where the first part hits a mature miRNA sequence (or some other small RNA / query RNA category, we expect these to be the first part of the chimera), and the second part hits the target genome. Because the second (i.e target) part can be quite short at times and can map to multiple locations in the genome, I use ShortStack's `--mmap u` mode to assign multi mapping reads.

Steps in the script [mirfirst_chimeras.pl]:

1. Find the reads where the first part matches a miRNA and the second part is long enough to keep:
    1. Use bowtie2 to map reads to the smallRNA database (eg mature miRNA). The bowtie2 index with that prefix should already exist.
    1. Select only those reads where the first part matches the mature miRNA, and the remaining part is longer than a minimum threshold (default 18bp)
    1. Output the second part of these reads (with modified headers) to a temporary fasta file called btmpXXXXXXX (where bXXXXXXX is randomly generated)

1. Find which part of the second part of the original read matches with the target genome, and whether it is low-complexity
    1. Before we do this, we have to first find WHICH part of the second part matches the target genome.
       This is because shortstack uses bowtie (not bowtie2) for mapping reads, and bowtie does not have a local-alignment mode.
       If the second part of the read is 80 bp long, but only bases 30-60 map to the target genome, then shortstack will not be able to find the match.
       Therefore the next step is to use bowtie2 (which has a local-alignment mode), to map btmpXXXXXXX against the target genome,
       pull out the part that matches the target genome, and store the output in a new temp fasta file called b2tmpXXXXXXX
    1. The b2tmpXXXXXXX temp fasta file may have many "dust" / low complexity sequences (eg "AAAAAAAAA", "ATATATATATATATATAT", "AAACAACAACAAC" etc, so we run it through dustmasker, and only keep those reads that are above a certain min length (set globally using $min_length, default 18), once the low complexity part has been removed. The output of this step is stored in a third temporary fasta file called b3tmpXXXXXXX

1. Create a multifasta file from the uniqfa file, and run ShortStack
    1. ShortStack works on individual reads, not on uniqfa sequences, so it needs to be deduped which is done inside this script using `uniqfa_to_multifa.pl`
    2. This multifa file is the input for ShortStack. The ShortStack file output directory has the same postfix as given in `-n` or `--name`

To run this script, you need to have ShortStack installed. You also need csvtk and bedtools, so perhaps the simplest way is to use conda:
```
conda create  -n bucklab
conda install -n bucklab -c bioconda shortstack csvtk bedtools
conda activate bucklab

for a in sample1 sample2 sample3 sample4
do
  echo $a
  cd /project/$a
  mirfirst_chimeras.pl \
    -i $a.cutadapt.uniq.fa \
    -s /databases/Mirbase/22/hsa-mature.fa \
    -t /databases/Gencode/H.sapiens/25/GRCh38.p7.genome.fa \
    -n $a.hsa-mature.human.18 \
    -m 18
  cd ..
done
```
Once the ShortStacks have finished running, we need to convert them to bam files, where the chimera read name is modified to include the miRNA name that the first part hit:
```
for a in sample1 sample2 sample3
do
  echo $a
  cd /project/$a

  bowtie2 --threads 8 --no-unal --no-head --norc -L 10 -i C,1 -N 1 --score-min C,30 --sensitive-local \
    -x /databases/Mirbase/22/hsa-mature.fa -f $a.cutadapt.uniq.fa \
  | extract_alignment_coordinates_from_sam.pl \
  > $a.cutadapt.uniq.fa.hsa-mature.tsv

  samtools view -F 4 ShortStack_noranmax.$a.hsa-mature.human.18.mall_mincov1/$a.hsa-mature.human.18.multi.bam | sed 's/:/\t/' \
  | csvtk join -T -t -H - <(cut -f1,5 $a.cutadapt.uniq.fa.hsa-mature.tsv) \
  | perl -plne 's/\t/:/; s/(.*)\t(\S+)$/$2:$1/' \
  | cat <(samtools view -H ShortStack_noranmax.$a.hsa-mature.human.18.mall_mincov1/$a.hsa-mature.human.18.multi.bam) - \
  | samtools view -b \
  > $a.hsa-mature.human.18.mirfirst.bam
done
```

## Get annotation files ready

Once we identify the 
