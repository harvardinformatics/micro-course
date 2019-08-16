Introduction to Bioinformatics
=========================

Much of bioinformatics these days involves working with high throughput sequencing data, typically in a UNIX environment. This workshop will give a quick refresher on the UNIX command line, introduce a number of common file formats and how to process them, and along the way give you some more advanced tips and tricks. The goal here is not to exhaustively cover each tool, but rather to show you some ways of solving common problems and along the way hopefully learn some useful shell commands.

For file formats, we’ll focus on:
 * fastq/fasta (the most common way to store sequence data, and what you’ll get back from the core after doing sequencing)
 * gff (the most common format for storing gene models and other complex genome annotations)
 * bed (a flexible format for storing any kind of interval data)
 * sam/bam (the dominant format for storing reference-based alignments of sequence data to a genome)

While some of these (e.g. fasta) are hopefully familiar to you, others may not be. You can find extensive detail about each file format from the following links:
 * [fastq](https://en.wikipedia.org/wiki/FASTQ_format)
 * [gff](https://useast.ensembl.org/info/website/upload/gff3.html)
 * [bed](https://genome.ucsc.edu/FAQ/FAQformat.html#format1)
 * [sam/bam](http://www.htslib.org/doc/sam.html)

In the morning, we'll cover fastq/fasta, gff, and bed files, as well as some useful Unix tricks. In the afternoon, we'll cover read mapping and sam/bam files.

Setup - Navigating the File System
----------

You'll need to copy some data to work with the examples below. We'll use this as an opportunity to quickly review navigating the file system using the command line.

First we want to make a new directory, using the `mkdir` command, to put the data you copy

```
mkdir -p intro_bioinf_2019
```

A few things to note: first, a general command in UNIX will be a program (`mkdir`) followed by options (`-p`). Options will either start with `-` or `--`, and many programs also have *positional* arguments, e.g. the directory we want to create in this case. We'll see this syntax repeatedly throughout the day.

Typing a program with the `-h` or `--help` option and nothing else will usually give you some basic info about it, e.g.:

```
mkdir --help
```

Now let's copy some data. We want to copy all the files from `/n/holylfs/LABS/informatics/workshops/bionano-gtt-data/` into a new subdirectory (`data`) in the directory we just made.

First we need to make the data subdirectory. We can do this in two ways:

Option 1
```
cd intro_bioinf_2019
mkdir data
```

Option 2
```
mkdir intro_bioinf_2019/data
```

In the first case, we first move into the directory with `cd` and then use mkdir to make a new directory in the current location.

In the second case, we don't move where we are on the file system, instead we specify path we want to create. Note that this is still *relative* to our current directory. We could do this instead with an *absolute* path, which starts from / (the root of the fileystem). We'll see this in the copy command. **NOTE** you may need to change this command depending on where in the filesystem you are. If you are not sure, use `pwd` to **p** rint **w** orking **d** irectory.

```
cp -v /n/holylfs/LABS/informatics/workshops/bionano-gtt-data/* intro_bioinf_2019/data
```

Now, let's change directory (using `cd`) so that we are in the `intro_bioinf_2019` directory (**not** the `data` subdir). Note that just typing `cd` with no arguments will bring you to your home directory, if you are confused about where you are in the file system and need to reset.

```
cd intro_bioinf_2019
```

Sequencing data (FASTQ/FASTA)
--------

To start with, let’s focus on a couple of things you might want to do with sequence data. One common task with fastq data is subsampling: you may want to do some preliminary analysis or tests with just a small subset of your sequencing reads in order to check to make sure things work before using your whole dataset. One way to do this is just take the first X sequences (which is easy to do with the head command), but best practice is to sample randomly from your file. We can do this with seqtk, which has a sample command.

Let's do this both ways. But first, we need to talk about output redirection. By default commands will often just print the output to the screen (stdout); instead we want to save this to a file, which we do using `>`, e.g.:

```
ls data > bioinf2019_data_files
```

Okay now back to getting the first 1000 lines. We'll see what the head command does by default.

```
head data/Falb_COL2.1.fastq
```

Each fastq record is 4 lines long, so to get 1000 of them we need 4000 lines. Using the help for head, how do you get head to return first 4000 lines of a file? Remember to use output redirection to save this to a file intead of print to your screen.

```
head -n 4000 data/Falb_COL2.1.fastq > Falb_COL2.1.subsample_head.fastq
```

Now let's use seqtk to actually sample

Get the help:
```
seqtk sample
```
Run the command:
```
seqtk sample -s 42 data/Falb_COL2.1.fastq 1000 > Falb_COL2.1.subsample_seqtk.fastq
```

For paired-end data, if we use the same seed for both files, we'll get the same output and our reads will remain paired:
```
seqtk sample -s 42 data/Falb_COL2.2.fastq 1000 > Falb_COL2.2.subsample_seqtk.fastq
```

Okay now we want to verify that we actually got 1000 records. For fastq files we can use the fact that a record is (almost always) four lines, and use the `wc` command. If you haven't used this command before, start with `wc --help`.

```
wc -l *.fastq
```

Another common operation on fasta files is extracting just a subset of a larger files, usually by name. For example, we may want to extract just a chromosome from a whole genome fasta file. There are a bunch of ways to do this; I’ll show you two.

Option 1: using seqtk

To use seqtk, we need to make a file with the name of the sequences we want to extract. This is inefficient for only a few regions, but very useful if we want to get a large number of sequences.

echo "X" > seqtk.regions
seqtk subseq data/dmel-all-chromosome-r6.20.fasta seqtk.regions > dmel-X.seqtk.fa
Option 2: using samtools

We can also use samtools, which is faster for just one region but inefficient for many sequences.

samtools faidx data/dmel-all-chromosome-r6.20.fasta X > dmel-X.samtools.fa

Now say we want to get the length of the Y chromosome in the Drosophila assembly; we are going to do this by stringing together three commands with pipes (`|`).

```
#extract Y sequence using samtools
```

```
#remove header with tail
```

```
count characters
```

```
combine with pipe
```

Interval and annotations data (BED, GFF)
--------

TO DO: INtroduce file types

<AWK COMES HERE>

Convert GFF to BED

Now, let’s say we want to extract just the X chromosome sequence and merge all the overlapping intervals in this CDS file, to give us a set of non-redundant intervals that contain all the protein-coding sequence on the X chromsome. We can do this using bedtools. Start by just typing bedtools to see all the subcommands.

In this case, merge looks like what we want. Note that many bedtools commands can take GFF as input, and will do the conversion for you automatically to bed:

grep "^X" hg38_cds.gff | bedtools merge -i - > hg38_chrX_cds.bed
Note that we get an error here, saying we need sorted input. Let’s try adding a sort step before the bedtools merge, using bedtools sort:

grep "^X" hg38_cds.gff | bedtools sort -i - | bedtools merge -i - > hg38_chrX_cds.bed
Woo hoo! It works.

Here is an exercise for you: let’s say we want a bed file that is the exact inverse of our hg38_chrX_cds.bed, that is all the regions of chromosome X that are not protein-coding (including UTRs, introns, etc). Look at the bedtools menu and figure out which tool could be used for this. Hint: you’ll need the hg38.genome file in your data directory, and you’ll also need to include a grep step (see if you can figure out why).

Another thing we might want to do is extract sequence data for particular intervals. You can do this in a couple of ways. By way of example we’ll use the dmel-genes.bed file. First, though, let’s pick 10 random genes to work with using bedtools sample, which selects random intervals from a file:

bedtools sample -i data/dmel-genes.bed -n 10 > dmel-10genes.bed
To get sequence for these genes, we can use bedtools getfasta:

bedtools getfasta -fi data/dmel-all-chromosome-r6.20.fasta -bed dmel-10genes.bed -name > dmel-10genes-bedtools.fasta
We need to use the name option to include the bed name field in the fasta definition line.

We can also extract sequence with seqtk, which can take a bed file as an interval list:

seqtk subseq data/dmel-all-chromosome-r6.20.fasta dmel-10genes.bed > dmel-10genes-seqtk.fa
Anything different between the files?

We are going to look at a few more uses of bedtools, but first let’s do an exercise that will show you how you can string a lot of commands together to accomplish a task. The goal is to get the DNA sequence for 10 random exons from chromosome 2L in D. melanogaster. Take this in steps, and you can do testing using temporary files but in the end you’ll need to link everything in one command.

e are going to move on to genome arithmetic now, that is doing calculations on intervals. One common question is, given a list of positions (e.g., peak calls or enhancers), what elements from another list (e.g. genes) are closest to each one. We’ll be using example files consisting of human-specific and shared with chimp neural crest enhancers published a few years ago. We can use bedtools closest for this:

bedtools closest -a data/human_enhancers.bed -b data/ucscGenes.bed | head
In this case, does bedtools output a bed file? How could we convert this to a proper bed file, where we keep the score but replace the names (HSE####) with the closest gene? (Hint: think about UNIX tools – e.g. cut).

Now let’s make a file, that has for each enhancer, the distance to the closest gene in the score file and the name of the closest gene in the name file. We’ll use the -t “first” option as well to guarantee only one line for each enhancer:

bedtools closest -a data/human_enhancers.bed -b data/ucscGenes.bed -d -t "first" | cut -f1,2,3,9,10 > hse_closest_genes.bed
Another common task is finding overlaps between two sets of intervals. For example, we might want to know which of our enhancers identified in neural crest cells are shared with enhancers identified in some other cell type. In your data directory you should see another bed file, H1-hESC-H3K27Ac.bed, which has H3K27Ac peaks from human embryonic stem cells. We’ll use bedtools to get the overlaps between this file and each of our enhancer files. The output of this can be a bit complicated, so let’s start by just looking at a few things:

bedtools intersect -a data/human_enhancers.bed -b data/H1-hESC-H3K27Ac.bed | head -n 3
What are we seeing here?

bedtools intersect -a data/human_enhancers.bed -b data/H1-hESC-H3K27Ac.bed -wa | head -n 3
bedtools intersect -a data/human_enhancers.bed -b data/H1-hESC-H3K27Ac.bed -wb | head -n 3
bedtools intersect -a data/human_enhancers.bed -b data/H1-hESC-H3K27Ac.bed -loj | head -n 3
bedtools intersect -a data/human_enhancers.bed -b data/H1-hESC-H3K27Ac.bed -wao | head -n 3
Note that all these options (-wa, -wb, -loj, -wao) change what information is returned, and can be very useful for different purposes.

Let’s do an exercise now. What if we wanted to get a new bed file with only the neural crest specific enhancers? That is, we want to filter out anything that overlaps (by even 1 bp) with a H3K27Ac peak in stem cells. Look at the bedtools intersect -h output and see if you can figure out what option to use to do this. Answer:

bedtools intersect -a data/human_enhancers.bed -b data/H1-hESC-H3K27Ac.bed -v > human_nconly_enhancers.bed
What if instead of removing the ones that overlapped, we just wanted to get data on how many stem cell enhancers each neural crest enhancer overaps?

Answer:

bedtools intersect -a data/human_enhancers.bed -b data/H1-hESC-H3K27Ac.bed -f 0.20 -r -c > human_overlaps.txt
Let’s say we wanted to shuffle enhancers to random locations in the genome, but exclude them from known genes and introns (just allow them to fall in intergenic regions. Biologically this is not completely realistic, since enhancers can be in introns, but it will illustrate the procedure. We might want to do this to, for example, ask whether enhancers are closer to genes than might be expected by chance. To do this, we need to do a few steps:

First, we make a bed file of intergenic regions (note we need to include a sort command):

bedtools sort -i data/ucscGenes.bed -faidx data/hg19.genome | bedtools complement -i - -g data/hg19.genome > hg19_intergenic.bed
Next, we use shuffle with the -incl option to tell bedtools where it can put our shuffled elements:

bedtools shuffle -incl hg19_intergenic.bed -i data/common_enhancers.bed \
-noOverlapping -g data/hg19.genome > common_enhancers_random.bed
We’ll end with a brief introduction to=some sam/bam processing.

One obvious task is to convert between sam and bam. We did this with Picard in the variant calling workshop, but we can also do it with samtools:

samtools view data/Falb_COL2.final.bam | head
samtools view -b data/Falb_COL2.final.sam > Falb_COL2.new.bam
We can also convert a bam file to a bed file, which can be useful for a number of reasons. We’ll do this in two steps: first, we need to produce a bam file that is sorted by query name (not start position), and excludes reads that are not properly paired:

samtools view -b -f 0x2 data/Falb_COL2.final.bam | samtools sort -n - | samtools fixmate -r - Falb_COL2.querysort.bam
Then we can use bamtobed to produce a bed file, where each interval is a read:

bedtools bamtobed -i Falb_COL2.querysort.bam -bedpe > Falb_COL2.bed
Finally, let’s look at calulating depth of coverage:

samtools depth data/Falb_COL2.final.bam > Falb_COL2.depth
bedtools genomecov -ibam data/Falb_COL2.final.bam -bga -pc -g data/Falb_region.genome > Falb_COL2_depth.bedGraph
