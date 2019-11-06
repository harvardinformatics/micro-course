Introduction to Bioinformatics
=========================

Much of bioinformatics these days involves working with high throughput sequencing data, typically in a UNIX environment. This workshop will give a quick refresher on the UNIX command line, introduce a number of common file formats and how to process them, and along the way give you some more advanced tips and tricks. The goal here is not to exhaustively cover each tool, but rather to show you some ways of solving common problems and along the way hopefully learn some useful shell commands.

We’ll focus on the following the following file formats (click on the links for extensive details about each file format):
 * [FASTQ](https://en.wikipedia.org/wiki/FASTQ_format)/[FASTA](https://en.wikipedia.org/wiki/FASTA_format) (the most common ways to store sequence data, and what you’ll get back from the core after doing sequencing)
 * [GFF (GFF3)](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md) (the most common format for storing gene models and other complex genome annotations)
 * [BED](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) (a flexible format for storing any kind of interval data)
 * [SAM/BAM](https://samtools.github.io/hts-specs/SAMv1.pdf) (the dominant formats for storing reference-based alignments of sequence data to a genome)

In this session, we'll cover fastq/fasta, gff, and bed files, as well as some useful Unix tricks.

Setup - Navigating the File System
----------

You'll need to copy some data to work with the examples below. We'll use this as an opportunity to quickly review navigating the file system using the command line.  

At the highest level in the structure is the root, signified by a forward slash `/`, which contains all lower directories (aka 'folders'). The path to all other directories 'lower' in the tree are also separated by a `/`. For example, to specify the path for a directory called 'bin' that is contained within the directory called 'usr' (i.e. two levels 'down' from root), we'd type `/usr/bin`. This is referred to as the **absolute path**, where the location of a directory/file is specified from the root directory.

If instead we didn't put the first `/` in the path, the shell would look in the *current directory* for a directory called 'usr', and would probably say that the directory doesn't exist. This is referred to as the **relative path**, where the path starts in the present working directory. Be careful not to get the two mixed up!

-----
There are several vital commands for navigating your way around a file system:
```
pwd: 'present working directory'. This prints where in your file system you currently are. If you ever get lost, use this!

ls: 'list'. This prints the contents of the directory you are in

cd: 'change directory'. This moves you to a different directory that you specify.

mkdir: 'make directory'. This creates a directory in your current directory or wherever you specify.
```

Now let's actually copy the data we need! First we want to make a new directory, using the `mkdir` command, to put the data you copy:

```
mkdir -p intro_bioinf_2019
```

A few things to note: first, since we don't specify the absolute path, it creates the directory within our current directory. Second, a general command in UNIX will be a program (`mkdir`) followed by options (`-p`, which means "create this directory/directories, including any parent directories, and don't error if any already exist").
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

Now let's copy the data. We want to copy all the files from `/n/holylfs/LABS/informatics/workshops/bionano-gtt-data/` into a new subdirectory (`data`) in the directory we just made.

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

In the second case, we don't move where we are on the file system, instead we specify path we want to create. Note that this is still *relative* to our current directory. We could do this instead with an *absolute* path, which starts from the root of the fileystem. We'll see this in the copy command. **NOTE** you may need to change this command depending on where in the filesystem you are. If you are not sure, use `pwd` to **p** rint **w** orking **d** irectory.

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

We can also use a different program called samtools, which is faster for just one region but inefficient for many sequences.

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


Manipulating files with grep
--------
We've explored a little bit with how to fish information out of sequence files using programs like seqtk and samtools, but we'd like something a little more robust. **grep** is a powerful command-line search tools that is included as part of most Unix-like systems. It is one of the most useful tools in bioinformatics! At the most basic level, grep searches for a string of characters that match a pattern and will print lines containing a match. Basic syntax is:  

`grep 'pattern' file_to_search`  


By default, grep will match any part of the string, so for example:

```
echo 'my dog is brown' > sample.txt  
grep 'dog' sample.txt  
grep 'do' sample.txt  
grep 'd' sample.txt
```


The above three grep commands will all match the text in the example file and will print the line. However, there are a huge number of arguments that can modify how grep behaves. Here are a few useful examples!


`grep -w` matches *entire words*.
- So in the above example:  
`grep -w 'dog' sample.txt` would match the string, but `grep -w 'do' sample.txt` would not.

`grep -i` allows case-insensitive matches
- In the above example, `grep -i 'DOG'` would still match the line

`grep -v` *inverts*, returning lines that *do not* match the pattern.

`grep -o` returns only the matching words, not the entire line.

`grep -c` counts the number of lines that match the pattern.
- Equivalent to `grep 'pattern' file | wc -l`

Print lines before/after a match:
- `grep -A [n]` returns matching line and *n* lines after match
- `grep -B [n]` returns matching line and *n* lines before match
- `grep -C [n]` returns matching line and *n* lines before *and* after match

Note: `grep -A` is very useful for pulling out specific lines from a FASTA file...just make sure your FASTA file is single line and not multi-line!

`grep -f pattern.txt` matches a list of patterns contained in pattern.txt against target file.
- I.e. the command `grep -f patterns.txt file.txt` will match all patterns in patterns.txt against file.txt
- Pattern file has one pattern per line


Remember that we can combine different grep arguments with each other! E.g.:

`grep -w -A 1 '>X' data/dmel-all-chromosome-r6.20.fasta`  
What will this line do??

`grep -w -v -f first3.txt data/human_enhancers.bed | less`  
How about this line?

There are many other functions of grep! When in doubt, remember you can check the help page using `man grep`...or just by using google :)

### Pattern matching with regular expressions
Regular expressions, aka "regex", are patterns that describe sets of strings. In other words, they allow you to match complex patterns with grep, not just exact matches. Regex is extremely powerful, but can also get (very) complicated, so we'll just stick to a few basic uses.

Regex has certain meta-characters that are reserved for special uses:
```
^: matches pattern at start of string
$: matches pattern at end of string
.: matches any character except new lines
[]: matches any of enclose characters
[^]: matches any characters *except* ones enclosed (note: is different from ^)
\: "escapes" meta-characters, allows literal matching
```

Note that if we want to match any of these special characters literally (e.g. matching a period character "."), we would need to use a "\\" to escape it first:  
`grep '\.' data/dmel-all-no-analysis-r6.20.gff` will literally match a period.

One example of how regex can come in handy is using the ^ special character to quickly count how many sequences are in a FASTA file, which we would do as follows:

`grep -c '^>' data/mel-all-chromosome-r6.20.fasta`

This command matches lines in the FASTA file that start with a ">" character, i.e. the header lines, and uses the -c argument to count how many matches!


### Practice:

- How many sequences are in the file Falb_COL2.1.fastq?

`grep -c '^@' data/Falb_COL2.1.fastq`

- Filter out lines matching unassigned contigs (chrUn) in the file hg19.genome and direct the output to a file.

`grep -v 'chrUn' data/hg19.genome  > hg19_noUn.genome`

- Use grep to pull out the header line of **only** the 2R chromosome arm from dmel-all-chromosome-r6.20.fasta and direct the output to a file.

`grep -w '^>2R' data/dmel-all-chromosome-r6.20.fasta > 2R_header.txt`

- Use grep to extract the lines of only the *major chromosome arms* (2L, 2R, 3L, 3R, and X) from the file dmel-all-no-analysis-r6.20.gff and pipe to less

```
printf '^>%s\n' 2L 2R 3L 3R X > major_arms.txt
grep -w -f major_arms.txt \
data/dmel-all-no-analysis-r6.20.gff | less
```
