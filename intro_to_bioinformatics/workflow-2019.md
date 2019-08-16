Introduction to Bioinformatics
=========================

Much of bioinformatics these days involves working with high throughput sequencing data, typically in a UNIX environment. This workshop will give a quick refresher on the UNIX command line, introduce a number of common file formats and how to process them, and along the way give you some more advanced tips and tricks. The goal here is not to exhaustively cover each tool, but rather to show you some ways of solving common problems and along the way hopefully learn some useful shell commands.

For file formats, we’ll focus on:
 * fastq/fasta (the most common way to store sequence data, and what you’ll get back from the core after doing sequencing)
 * gff (the most common format for storing gene models and other complex genome annotations)
 * bed (a flexible format for storing any kind of interval data)
 * sam/bam (the dominant format for storing reference-based alignments of sequence data to a genome)

While some of these (e.g. fasta) are hopefully familiar to you, others may not be. You can find extensive detail about each file format from the following links:
 * fastq[https://en.wikipedia.org/wiki/FASTQ_format]
 * gff[https://useast.ensembl.org/info/website/upload/gff3.html]
 * bed[https://genome.ucsc.edu/FAQ/FAQformat.html#format1]
 * sam/bam[http://www.htslib.org/doc/sam.html]

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

```
cd intro_bioinf_2019
mkdir data
```

```
mkdir intro_bioinf_2019/data
```

In the first case, we first move into the directory with `cd` and then use mkdir to make a new directory in the current location.

In the second case, we don't move where we are on the file system, instead we specify path we want to create. Note that this is still *relative* to our current directory. We could do this instead with an *absolute* path, which starts from / (the root of the fileystem). We'll see this in the copy command. *NOTE* you may need to change this command depending on where in the filesystem you are. If you are not sure, use `pwd` to *p* rint *w* orking *d* irectory.

```
cp -v /n/holylfs/LABS/informatics/workshops/bionano-gtt-data/* intro_bioinf_2019/data
```

Sequencing data (FASTQ/FASTA)
--------


Gene models and genome annotations (GFF)
--------


Interval data (BED)
--------


Alignments to a Reference (BAM/SAM)
----------------
