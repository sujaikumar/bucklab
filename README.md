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

Get the bowtie2 databases (for quickly searching) and bowtie database (for shortstack) ready. hsa-mature.fa should be a DNA file, not RNA:
```
bowtie2-build /databases/Mirbase/22/hsa-mature.fa /databases/Mirbase/22/hsa-mature.fa
bowtie2-build /databases/Gencode/H.sapiens/25/GRCh38.p7.genome.fa /databases/Gencode/H.sapiens/25/GRCh38.p7.genome.fa
bowtie-build  /databases/Gencode/H.sapiens/25/GRCh38.p7.genome.fa /databases/Gencode/H.sapiens/25/GRCh38.p7.genome.fa
```
Steps in the script [mirfirst_chimeras.pl](mirfirst_chimeras.pl)

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
    1. ShortStack works on individual reads, not on uniqfa sequences, so it needs to be deduped which is done inside this script using [uniqfa_to_multifa.pl](uniqfa_to_multifa.pl)
    2. This multifa file is the input for ShortStack. The ShortStack file output directory has the same postfix as given in `-n` or `--name`

To run this script, you need to have ShortStack installed. You also need csvtk and bedtools, so perhaps the simplest way is to use conda:
```
conda create -n bucklab -c bioconda shortstack csvtk bedtools bowtie2 bowtie cutadapt blast parallel
conda activate  bucklab

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
for a in sample1 sample2 sample3 sample4
do
  echo $a
  cd /project/$a

  bowtie2 --threads 8 --no-unal --no-head --norc -L 10 -i C,1 -N 1 --score-min C,30 --sensitive-local \
    -x /databases/Mirbase/22/hsa-mature.fa -f $a.cutadapt.uniq.fa \
  | extract_alignment_coordinates_from_sam.stranded.pl \
  > $a.cutadapt.uniq.fa.hsa-mature.tsv

  samtools view -F 4 ShortStack_noranmax.$a.hsa-mature.human.18.mall_mincov1/$a.hsa-mature.human.18.multi.bam | sed 's/:/\t/' \
  | csvtk join -T -t -H - <(cut -f1,5 $a.cutadapt.uniq.fa.hsa-mature.tsv) \
  | perl -plne 's/\t/:/; s/(.*)\t(\S+)$/$2:$1/' \
  | cat <(samtools view -H ShortStack_noranmax.$a.hsa-mature.human.18.mall_mincov1/$a.hsa-mature.human.18.multi.bam) - \
  | samtools view -b \
  > $a.hsa-mature.human.18.mirfirst.bam
done
```
What this is doing is prepending each read in the bam file with the miRNA name that the first part hits

Next Steps:

## Convert bam files to miRNA-target-loci csv tables with counts from each sample

First we need a merged bed file of all the miRNA-specific target loci where the chromosome name has the miRNA name attached:
```
# this $OUTPREFIX prefix will be used for all subsequent output files. I'm just using test for now
OUTPREFIX=test

for a in sample1 sample2 sample3 sample4
do
  samtools view $a/hsa-mature.human.18.mirfirst.bam \
  | extract_alignment_coordinates_from_sam.stranded.pl \
  | perl -lne '
    if (/^(\S+?):\S+\t\d+\t\d+\t\d+\t(\S+)\t(\d+)\t(\d+)\t\d+\t(.)$/) {
      ($mirfirst, $chr, $st, $en, $strand) = ($1, $2, $3, $4, $5);
      print "$mirfirst:$chr\t$st\t$en\t.\t.\t$strand"
    }'
done \
| sort -k1,1 -k2,2n \
| bedtools merge -s -i - -c 6 -o distinct | perl -lne 'print "$1\t$2\t$3\t$1:$2-$3:$4\t.\t$4" if /^(\S+)\t(\d+)\t(\d+)\t(.)/' \
> $OUTPREFIX.bed
```

Now, we go back to the bam files and create a count for each sample for each of these bed regions above:
```
for a in sample1 sample2 sample3 sample4
do
  samtools view $a/$a.hsa-mature.human.18.mirfirst.bam \
  | extract_alignment_coordinates_from_sam.stranded.pl \
  | perl -lne 'print "$2:$3\t$4\t$5\t$1\t.\t$6" if /^((\S+?):\S+)\t\d+\t\d+\t\d+\t(\S+)\t(\d+)\t(\d+)\t\d+\t(.)$/' \
  | intersectBed -s -a - -b ../$OUTPREFIX.bed -bed -wa -wb \
  | awk -v sample=$a -F"\t" '{if(h[$4]++<1){print $10"\t"sample}}'
done \
| perl -ane '
    chomp;
    @F=split/\t/;
    $h{$F[0]}{$F[1]}++;
    $samples{$F[1]}++;
  }
  {
    for $mirfirstcluster (keys %h) {
      print $mirfirstcluster;
      for $sample (sort keys %samples) {
        print "\t" . $h{$mirfirstcluster}{$sample};
      }
      print "\n"
    }
' | pigz -c > $OUTPREFIX.samplecounts.tsv.gz
```

## Get annotation files ready

I use https://www.gencodegenes.org for human and mouse because they have well formatted files. Broadly - the steps are to download the relevant files from Gencode, and then to create separate files for each annotation type (UTR3, CDS, UTR5, rRNA, tRNA, miRNA). I also use the ensembl regulatory regions files and get the separate types of reg regions in separate files

```
# am using a bash variable here so that new versions just require a one variable change
r=34
s=human

# or you could do a mouse one:
# r=M25
# s=mouse

mkdir -p Gencode/$s/$r
cd !$

wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_$s/release_$r/gencode.v$r.chr_patch_hapl_scaff.annotation.gff3.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_$s/release_$r/gencode.v$r.tRNAs.gff3.gz

g=gencode.v$r.chr_patch_hapl_scaff.annotation.gff3.gz
t=gencode.v$r.tRNAs.gff3.gz

mkdir genomespace
cd    genomespace
zcat                                 ../$t | perl -plne 's/gene_n.*transcript_name=([^;]+)/gene_name=$1/' > tRNA.gff3
zgrep -P "\texon\t.*type=miRNA"      ../$g | perl -lne '@F=split/\t/;$F[2]="miRNA";  print join("\t",@F)' > miRNA.gff3
zgrep -P "\texon\t.*type=(Mt_)?rRNA" ../$g | perl -lne '@F=split/\t/;$F[2]="rRNA";   print join("\t",@F)' > rRNA.gff3
zgrep -P "\texon\t.*type=lincRNA"    ../$g | perl -lne '@F=split/\t/;$F[2]="lncRNA"; print join("\t",@F)' > lncRNA.gff3

zgrep -P "\tfive_prime_UTR\t"  ../$g    > UTR5.gff3
zgrep -P "\tthree_prime_UTR\t" ../$g    > UTR3.gff3
zgrep -P "\tCDS\t"             ../$g    > CDS.gff3

# the original Gencode GFF files don't have introns,
# and I found this tool/method from the agat package
# as the best/easiest way for getting (protein-coding) introns from the gff3:

conda install -c bioconda agat 

zgrep -P "type=protein_coding" ../$g    > pc.gff3
agat_sp_add_introns.pl --gff pc.gff3 --out introns.gff
grep -P "\tintron\t" introns.gff | sort -k1,1 -k4,4n > intron.gff3
```
## Regulatory regions
Download regulatory regions from ensembl (make sure they match the genome version, GRCh38 matches GencodeGenes human releases 20-34, and Mouse GRCm38 matches mouse releases M2-M25 - see https://www.gencodegenes.org/human/releases.html):

```
cd Gencode/$s/$r
wget ftp://ftp.ensembl.org/pub/release-100/regulation/homo_sapiens/homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20190329.gff.gz
wget ftp://ftp.ensembl.org/pub/release-100/regulation/mus_musculus/mus_musculus.GRCm38.Regulatory_Build.regulatory_features.20180516.gff.gz
```
The mouse one doesn't have column 3 equal to the type of reg region it is, so that needs to be fixed first:
```
zcat mus_musculus.GRCm38.Regulatory_Build.regulatory_features.20180516.gff.gz \
| perl -lne '@F=split/\t/;if($F[8]=~/feature_type=([^;]+)/){$F[2] = $1; $F[2]=~s/ /_/g; print join("\t",@F)}' \
| gzip -c > mus_musculus.GRCm38.Regulatory_Build.regulatory_features.20180516.fixed.gff.gz
```

Now we can create sep annotation files for each type of reg region in the genomespace folder:
```
reg=mus_musculus.GRCm38.Regulatory_Build.regulatory_features.20180516.fixed.gff.gz
#reg=homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20190329.gff.gz

cd Gencode/$s/$r/genomespace

zcat ../$reg | cut -f3 | sort | uniq \
| parallel "zcat ../$reg | awk '\$3==\"{}\"' | perl -plne 's/^([1-9MX])/chr\$1/' > {}.gff3"
```

## Annotate each target site with protein-coding and noncoding features and regulatory regions

intersect all the mirfirst-target loci / clusters with all the annotation gff3 from previous step:

Note: These have to be done separately because reg regions don't have strand, so the bedtools intersect command is different for protein-coding vs reg regions

Note: For mouse, the names of the regulatory region .gff3 s are slightly different than for human

```
# cd back to directory with $OUTPREFIX.bed
# the miRNA-target loci in $OUTPREFIX.bed have the miRNA/smallRNA name attached to the chromosome name,
# so we have to first separate that to make it a proper BED file using perl -plne 's/.+?://'

parallel "
  bedtools intersect -s -wa -wb \
    -a <(perl -plne 's/.+?://' $OUTPREFIX.bed) \
    -b Gencode/$s/$r/genomespace/{}.gff3 \
  | cut -f4,15 | perl -plne 's/\\S*gene_name=([^;]+).*/\$1/' | awk 'a[\$1]++<1' \
  > $OUTPREFIX.{}.tsv" \
::: UTR3 CDS UTR5 miRNA intron tRNA rRNA lncRNA

parallel "
  bedtools intersect -wa -wb \
    -a <(perl -plne 's/.+?://' $OUTPREFIX.bed) \
    -b Gencode/$s/$r/genomespace/{}.gff3 \
  | cut -f4,15 | perl -plne 's/ID=[^:]+:([^;]+).*/\$1/' \
  | awk 'a[\$1]++<1' \
  > $OUTPREFIX.{}.tsv" \
::: CTCF_binding_site TF_binding_site enhancer open_chromatin_region promoter_flanking_region promoter

# for mouse regulatory regions it will be:
# ::: CTCF_Binding_Site TF_binding_site Enhancer Open_chromatin Promoter_Flanking_Region Promoter
```

## Find seed regions

## Run RNAhybrid between the query and target parts of each miRNA-targetsite pair
