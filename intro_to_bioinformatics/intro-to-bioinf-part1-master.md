Introduction to Bioinformatics
=========================

Much of bioinformatics these days involves working with high throughput sequencing data, typically in a UNIX environment. This workshop will give a quick refresher on the UNIX command line, introduce a number of common file formats and how to process them, and along the way give you some more advanced tips and tricks. The goal here is not to exhaustively cover each tool, but rather to show you some ways of solving common problems and along the way hopefully learn some useful shell commands.

We’ll focus on the following the following file formats (click on the links for extensive details about each file format):
 * [FASTQ](https://en.wikipedia.org/wiki/FASTQ_format)/[FASTA](https://en.wikipedia.org/wiki/FASTA_format) (the most common ways to store sequence data, and what you’ll get back from the core after doing sequencing)
 * [GFF (GFF3)](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md) (the most common format for storing gene models and other complex genome annotations)
 * [BED](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) (a flexible format for storing any kind of interval data)
 * [SAM/BAM](https://samtools.github.io/hts-specs/SAMv1.pdf) (the dominant formats for storing reference-based alignments of sequence data to a genome)

In the morning, we'll cover fastq/fasta, gff, and bed files, as well as some useful Unix tricks. In the afternoon, we'll cover read mapping and sam/bam files.

Setup - Navigating the File System
----------

You'll need to copy some data to work with the examples below. We'll use this as an opportunity to quickly review navigating the file system using the command line.

First we want to make a new directory, using the `mkdir` command, to put the data you copy:

```
mkdir -p intro_bioinf_2019
```

A few things to note: first, a general command in UNIX will be a program (`mkdir`) followed by options (`-p`, which means "create this directory/directories, including any parent directories, and don't error if any already exist").
Options will either start with:
* `-` (for short/single-letter options, which can be combined; e.g., `mkdir -pv` is equivalent to `mkdir -p -v` ), or
* `--` (for long/multi-character options; e.g. `mkdir --parents --verbose`)

Many programs also have *positional* arguments (aka operands), e.g. the directory we want to create in this case.
We'll see this syntax repeatedly throughout the day.

Typing a program with the `-h` or `--help` option and nothing else will usually give you some basic info about it, e.g.:

```
mkdir --help
```

*NOTE: this will not work on macOS; use the [man pages](https://en.wikipedia.org/wiki/Man_page) instead (e.g., `man mkdir`)*

Now let's copy some data. We want to copy all the files from `/n/holylfs/LABS/informatics/workshops/bionano-gtt-data/` into a new subdirectory (`data`) in the directory we just made.

First we need to make the data subdirectory:

Option 1
```
cd intro_bioinf_2019
mkdir data
```

Option 2
```
mkdir intro_bioinf_2019/data
```

In the first case, we first move into the directory with `cd` and then use mkdir to make a new directory in the current working directory.

In the second case, we don't move where we are on the file system, instead we specify path we want to create. Note that this is still *relative* to our current directory. We could do this instead with an *absolute* path, which starts from / (the root of the fileystem). We'll see this in the copy command. **NOTE** you may need to change this command depending on where in the filesystem you are. If you are not sure, use `pwd` to **p** rint **w** orking **d** irectory.

Assuming you executed `mkdir intro_bioinf_2019/data`:
```
cp -v /n/holylfs/LABS/informatics/workshops/bionano-gtt-data/* intro_bioinf_2019/data
```

*The ` *` is a shell pattern that matches any filename, while the `-v` option causes `cp` to show which file is currently being copied (otherwise, `cp` is silent until it completes).*

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

Each fastq record is 4 lines long, so to get 1000 of them we need 4000 lines.

```
head -n 4000 data/Falb_COL2.1.fastq > Falb_COL2.1.subsample_head.fastq
```

To efficiently view (but not edit) this file, we can use the [less](https://en.wikipedia.org/wiki/Less_(Unix)) command:
```
less Falb_COL2.1.subsample_head.fastq
```
The arrow keys can be used to scroll forward and backward in the file, while `/` can be used to search for a pattern.
To exit, type `q`.

`less` does not load the entire file into memory before displaying, and can be used to scroll through extremely large files.

Now let's use seqtk to actually sample.

Get the help:
```
seqtk sample
```
> ```
> Usage:   seqtk sample [-2] [-s seed=11] <in.fa> <frac>|<number>
>
> Options: -s INT       RNG seed [11]
         -2           2-pass mode: twice as slow but with much reduced memory
> ```

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
> ```
>  4000 Falb_COL2.1.subsample_head.fastq
>  4000 Falb_COL2.1.subsample_seqtk.fastq
>  4000 Falb_COL2.2.subsample_seqtk.fastq
> 12000 total
> ```

Another common operation on fasta files is extracting just a subset of a larger files, usually by name. For example, we may want to extract just a chromosome from a whole genome fasta file. There are a bunch of ways to do this; I’ll show you two.

Option 1: using seqtk

To use seqtk, we need to make a file with the name of the sequences we want to extract. This is inefficient for only a few regions, but very useful if we want to get a large number of sequences.

```
echo "X" > seqtk.regions
seqtk subseq data/dmel-all-chromosome-r6.20.fasta seqtk.regions > dmel-X.seqtk.fa
```

Option 2: using samtools

We can also use samtools, which is faster for just one region but inefficient for many sequences.

```
samtools faidx data/dmel-all-chromosome-r6.20.fasta X > dmel-X.samtools.fa
```

Now say we want to get the length of the Y chromosome in the Drosophila assembly; we are going to do this by stringing together three commands with pipes, which feeds the *standard output* of the command on the left of the `|` to the *standard input* of the command on the right.

```
samtools faidx data/dmel-all-chromosome-r6.20.fasta Y > dmel-Y.samtools.fa
```

```
tail -n +2 dmel-Y.samtools.fa > dmel-Y-seqonly.fa
```

```
wc -c dmel-Y-seqonly.fa
```

```
samtools faidx data/dmel-all-chromosome-r6.20.fasta Y | tail -n +2 | wc -c
```

Pipes are a *very* useful optimization when working with large files, they avoid needing to create intermediate files, and the commands in a pipeline execute concurrently.
Many bioinformatics utilities can read their input from pipes.

Interval and annotations data (BED, GFF)
--------

Now, we are going to move on from fasta files (although we'll come back to them in a bit) and introduce files with data stored as genomic intervals. Of course this includes things like gene models, but also the output of peak callers, enhancer predictions, and many other things. We’ll focus on two file formats, gff and bed. In both cases intervals are represented as lines in the file with a chromosome, start position, and end position specified. Confusingly, bed files use a 0-based start (first base of a chromosome is 0), and a 1-based end, while GFF files use a consistent 1-based interval. This is clearer with an example:

The interval from the first base to the hundredth base of a chromosome would be represented as start = 0, end = 100 in a bed file, and start = 1, end = 100 in a GFF. Another way to think about this is that bed files are 0-based, half-open (the start is included in the interval but not the end), while GFF are 1-based inclusive (both start and end are included). This is super-confusing and will take a while to get used to.

There are a few implications of this. For example, to compute length for a bed interval is just end-start, while for a GFF it is end-start+1. Also, a one-base interval (e.g., the position of a variable site) is encoded as start=end in a GFF file, but as start=end-1 in a bed file (e.g., a SNP at position 57 would be 57,57 in GFF but 56,57 in bed). Unfortunately the other implication is that off-by-one errors are common especially when you are not used to thinking about this. Sorry, this is just life in bioinformatics.

Bed and GFF files also have different fields. A bed file is: seqname, start, end, name, score, strand, thickStart, thickEnd, itermRgb, blockCount, blockSizes, blockStarts (the last several are display attributes if you load a bed on a genome browser). Only the first three fields are required, but you cannot skip fields (so if you want to include strand, you need to also include name and score in that order). A gff file is: seqname, source, feature, start, end, score, strand, frame, attribute and all fields are required.

Let’s start by exploring a GFF file. Use `less` to look at the human gene annotation GFF stored as: data/Homo_sapiens.GRCh38.91.gff3. You can see there are a bunch of different feature fields. Let's say we wanted a list of just the feature fields. We can use the Unix command `cut` to extract a particular field from a file.

```
head -n 1000 data/Homo_sapiens.GRCh38.91.gff3 | cut -s -f 3 | less
```

What we really want is actually just the unique values. To do this, we can combine the Unix commands `sort` and `uniq`:

```
cut -s -f 3 data/Homo_sapiens.GRCh38.91.gff3 | sort | uniq -c
```

`uniq` collapses duplicate adjacect lines into a single copy, which is prefixed by a count via the `-c` option.
The lines are first `sort`ed so each count reflects the the number of occurrences of a line in the entire input.

Finally, let's briefly look at a BED file before moving on to some more advanced Unix tricks that will help us work with these files.

```
less data/dmel-genes.bed
```

Note the differences with a GFF. As an exercise, let's use the same command we used above, but this time to get all the chromosomes (first field, so `-f 1`) in the dmel-genes bed file.

```
cut -f 1 data/dmel-genes.bed | sort | uniq -c
```

This was a very quick introduction. We are going to do more with BED and GFF files, but first we are going to introduce some more complex Unix tools that will help us greatly in working with these kinds of files.
