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

Manipulating files with AWK
-----

Invented in the 1970's, [awk](https://en.wikipedia.org/wiki/AWK) is a scripting language included in most Unix-like operating systems. It specializes in one-liner programs and manipulating text files.

In many cases, if you're parsing information from a text file (such as a BED file, FASTA file, etc.), you could write a Python script...or you could do it with awk in a single line!

### Syntax
awk scripts are organized as:

`awk 'pattern { action; other action }' file`

Meaning that every time that the pattern is true, awk will execute the action in the brackets. If no pattern is specified, the action will be taken for every line in the input file, e.g. the following command prints every line:

`awk '{print}' data/hg38.genome | less`

The two most important patterns are `BEGIN` and `END`, which tell the action to take place before any lines are read and after the last line.

 `awk 'BEGIN{sum=0} {sum+=1} END {print sum}' data/hg38.genome`

 The above line sets a variable at the start of the script, adds 1 to it every line, then prints its value at the end.

If a variable hasn't been initialized, it is treated as 0 in numeric expressions, and an empty string in string expressions&mdash;awk will not print an error!
So the following awk script also prints the number of lines in the file data/hg38.genome:

`awk '{sum+=1} END {print sum}' data/hg38.genome`

### Input and output
Input to awk is split into **records** and **fields**.
- By default, **records** are separated by newline character, i.e # of records = # of lines in input file
- Each record is subdivided into **fields**, i.e. columns, as determined by the field separator (see below)

There are several important built-in variable in awk. The fields (columns) of each record are referred to by `$number`, so the first column would be `$1`, second would be `$2`, etc. `$0` refers to the entire record.<br/>
So to print the second column of each line in the file, we'd use:

`awk '{print $2}' data/hg38.genome | less`

And if we wanted to print the second then the first:

`awk '{print $2,$1}' data/hg38.genome | less`

Note that when the different fields are separated with commas in the `print` statement, they are joined by the output field separator (the **OFS** variable, described below), which is by default a space.
If the comma is omitted between fields (e.g., `awk '{print $2 $1}'`, they are concatenated without a separator.
<br/>
<br/>

We can also print strings using using quotation marks:

`awk '{print "First column:" $1}' data/hg38.genome | less`

Which for every line of the file will print the text "First column:" followed by the value in the first field.

---
awk has several other built-in variables that are very useful for parsing text, including:

|  |   |
---|---|
| **FS** | field separator (default: white space) |
| **OFS** |  output field separator, i.e. what character separates fields when printing|
| **RS** | record separator, i.e. what character records are split on (default: new line) |
| **ORS** | output record separator |
| **NR** | number of records in input (# lines by default) |

Using these, we can convert between file formats, e.g. make a comma-separated text file into a tab-separated file:

`awk 'BEGIN{FS="," ; OFS="\t"} {print $0}' data/enhancers.csv > data/enhancers.tsv`


### Conditionals and pattern matching
Like other programming languages, awk allows conditional matching with if/else statements.

`awk '{if(condition) action; else other_action}'`

awk uses the following conditional operators:

| | |
|-|-|
|==|equal to|
|!=|not equal to|
|>|greater than|
|>=|greater than or equal to|
|<|less than|
|<=|less than or equal to|
|&&|AND|
| \|\| |OR|
| ! | NOT |

In addition, awk also supports string matching using regular expressions, using the following expressions:

| | |
|-|-|
|\~|matches|
|!~|does not match|

For string matching, the pattern being matched must be enclosed by slashes, like so:

`awk '{if($1 ~ /pattern/) print}'`

Note that if an action isn't specified, the default action is `{print}`, so the previous awk command is equivalent to the following, which specifies only a pattern expression:

`awk '$1 ~ /pattern/'`

---
### Example uses:

- Count number of sequences in a FASTQ file:

`awk 'END{print NR/4}' data/Falb_COL2.1.fastq`

**Note**: this is technically safer than the command we looked at earlier using grep, as yo don't have to worry about accidentally counting the quality line.

- Only print annotations on a specific scaffold (chr1) that fall between 1Mb and 2Mb:

`awk 'BEGIN{FS="\t";OFS="\t"} {if($1 == "chr1" && $2 >=1000000 && $2 <= 2000000) print}' data/ucscGenes.bed | less`

**Note**: when we specify that we only want annotations from chr1, we're using exact match (`== "chr1"`) and not pattern match (`~ /chr1/`)...why is this??

- Only print lines of GFF file that match the string "exon" in their third column:

`awk 'BEGIN{FS="\t"} {if($3 ~ /exon/) print $0}' data/dmel-all-no-analysis-r6.20.gff | less`

- Convert from GFF (genome feature file) to BED file

`grep -v '^#' data/dmel-all-no-analysis-r6.20.gff | awk 'BEGIN{FS="\t"; OFS="\t"} {print $1,$4-1,$5}' | less`

**Note**: remember that BED and GFF files have different coordinate systems, i.e. BED start coordinate is 0 based, half-open, GFF is 1-based inclusive! Also, we are first using grep to skip the header lines in the GFF file.

Alternatively, you could do the whole thing with only awk, no grep required:

`awk 'BEGIN{FS="\t"; OFS="\t"} !/^#/ {print $1,$4-1,$5}' data/dmel-all-no-analysis-r6.20.gff | less`

With this command, `!/^#/` is a pattern (like `BEGIN` or `END`) that tell awk to execute the print statement when the start of the line does not match a `#`. Use whichever makes the most sense to you!


### Practice
Using awk:<br/>
- Pull out only the CDS annotations from the GFF file dmel-all-no-analysis-r6.20.gff and output them in BED format

`awk 'BEGIN{FS="\t"; OFS="\t"} {if($3 ~ /CDS/) print $1,$4-1,$5}' data/dmel-all-no-analysis-r6.20.gff`

- Extract FASTA information from the SAM file Falb_COL2.final.sam (hint: you will need to remove the header lines first, they start with `@`!) (Another hint: sequence ID = 1st column, sequence = 10th column)

`awk 'BEGIN{FS="\t"; OFS="\n"} !/^@/ {print ">"$1,$10}' data/Falb_COL2.final.sam | less`

- Write a command to calculate that average of the 5th column (i.e. mapping quality score) of a tab-separated SAM file Falb_COL2.final.sam and output it (hint: as above, need to remove headers)

`awk 'BEGIN{FS="\t"; sum=0} !/^@/ {sum+=$5} END{print sum/NR}' data/Falb_COL2.final.sam`

- Calculate the average length of gene annotations *only on the 2L arm* from the file dmel-genes.bed (hint: you'll need to use a combination of grep and awk...)

`awk 'BEGIN{FS="\t"; sum=0} $1 == "2L" {len=$3-$2; sum=sum+len} END{print sum/NR}' data/dmel-genes.bed`


Genomic Ranges
============

Interval files are sufficiently common to have a number of tools developed specifically for interacting with them. In the next hour, we will introduce some of the power of one tool, `bedtools`. This program has a huge number of options and can do a lot of stuff. We'll only briefly introduce the full power of interval manipulation today.

Finding the closest gene to something
-------------

One very common task in interval manipulation is trying to find the closest interval from one file to intervals in another file, e.g. you have a list of putative ehnancers (as a bed file) and you want to know the closet gene to each one.  **Exercise**: Scroll through the list of options and figure out what the best sub-command would be to do this.

To get help on a specific command, type `bedtools closest -h`

You'll see that for many commands, bedtools works by taking a 'query' bed file (specified with the -a option) and comparing it to a 'target' bed file (specified with the -b option). Although for some tools it doesn't matter which you put as -a and which as -b, in some cases it does matter so it is worth thinking through a bit.

**Exercise**: to find which gene is closest to each enhancer in a given BED file, what should be specify as -a and what as -b?

Answer:
```
bedtools closest -a data/human_enhancers.bed -b data/ucscGenes.bed | head

chr1	1015066	1015266	HSE897	9.5706324138	chr1	1017197	1051736	uc001acu.2
chr1	1590473	1590673	HSE853	17.3206898329	chr1	1571099	1655775	uc001agv.1
chr1	2120861	2121064	HSE86	66.0424471614	chr1	2115898	2126214	uc031pkt.1
chr1	6336418	6336624	HSE394	14.8892086906	chr1	6324331	6453826	uc001amt.3
chr1	7404594	7404794	HSE315	24.228051023	chr1	6845383	7829766	uc001aoi.3
chr1	11941325	11941525	HSE322	17.5192630328	chr1	11917520	11918992	uc001atj.3
chr1	15055555	15055755	HSE962	11.7065788965	chr1	14925212	15444544	uc001avm.4
chr1	15478024	15478224	HSE354	17.5771586196	chr1	15438310	15478960	uc009voh.3
chr1	16065307	16065508	HSE264	23.9092446094	chr1	16062808	16067884	uc001axb.1
chr1	23244249	23244449	HSE434	17.1331422162	chr1	23243782	23247347	uc001bgg.1
```
**Exercise**: in this case, does bedtools output a bed file? How could we convert this to a proper bed file, where we keep the score but replace the names (HSE####) with the closest gene? (Hint: think about UNIX tools)

Answer: `bedtools closest -a data/human_enhancers.bed -b data/ucscGenes.bed | awk -v OFS='\t' '{print $1, $2, $3, $9, $5}' | head`

Okay, now let's make a file that contains the closest gene for each of the human-specific enhancers, in proper bed format, but we'll add two additional options (refer back to bedtools closest -h for the full option list): -d to give us the distance to the nearest gene, which we'll put in the 'score' field of the output bed file, and -t "first" means that if there are ties, we'll just keep one (so that we have each enhancer assigned to a single gene).

`bedtools closest -a data/human_enhancers.bed -b data/ucscGenes.bed -d -t "first" | awk -v OFS='\t' '{print $1, $2, $3, $9, $10}' > hse_closest_genes.bed`

**Exercise**: What do we need to change to get the equivalent set but for the common enhancers?

Answer: The input file name for bed file 'A', the output file name, and the field numbers for the gene and score.
```
bedtools closest -a data/common_enhancers.bed -b data/ucscGenes.bed -d -t "first" | awk -v OFS='\t' '{print $1, $2, $3, $7, $8}' > ce_closest_genes.bed
```
Finding overlaps between intervals
-------------

Another common task is finding overlaps between two sets of intervals. For example, we might want to know which of our enhancers identified in neural crest cells are shared with enhancers identified in some other cell type. In your directory you should see another bed file, H1-hESC-H3K27Ac.bed, which has H3K27Ac peaks from human embryonic stem cells. We'll use bedtools to get the overlaps between this file and each of our enhancer files.

The tool we'll use is intersect. Let's start by looking at the options: ```bedtools intersect -h```

The output of this can be a bit complicated, so let's start by just looking at a few things:

```
bedtools intersect -a data/human_enhancers.bed -b data/H1-hESC-H3K27Ac.bed | head -n 3

chr1	1590473	1590673	HSE853	17.3206898329
chr1	6336418	6336624	HSE394	14.8892086906
chr1	7404594	7404794	HSE315	24.228051023
```

What are we seeing here?

```
bedtools intersect -a data/human_enhancers.bed -b data/H1-hESC-H3K27Ac.bed -wa | head -n 3

chr1	1590473	1590673	HSE853	17.3206898329
chr1	6336418	6336624	HSE394	14.8892086906
chr1	7404594	7404794	HSE315	24.228051023

bedtools intersect -a data/human_enhancers.bed -b data/H1-hESC-H3K27Ac.bed -wb | head -n 3
chr1	1590473	1590673	HSE853	17.3206898329	chr1	1590264	1590895	.	1000	.
chr1	6336418	6336624	HSE394	14.8892086906	chr1	6160044	6490440	.	291	.
chr1	7404594	7404794	HSE315	24.228051023	chr1	7400919	7412366	.	279	.
```

What about now?

```
bedtools intersect -a data/human_enhancers.bed -b data/H1-hESC-H3K27Ac.bed -loj | head -n 3
chr1	1015066	1015266	HSE897	9.5706324138	.	-1	-1	.	-1	.
chr1	1590473	1590673	HSE853	17.3206898329	chr1	1590264	1590895	.	1000	.
chr1	2120861	2121064	HSE86	66.0424471614	.	-1	-1	.	-1	.

bedtools intersect -a data/human_enhancers.bed -b data/H1-hESC-H3K27Ac.bed -wao | head -n 3
chr1	1015066	1015266	HSE897	9.5706324138	.	-1	-1	.	-1	.	0
chr1	1590473	1590673	HSE853	17.3206898329	chr1	1590264	1590895	.	1000	.	200
chr1	2120861	2121064	HSE86	66.0424471614	.	-1	-1	.	-1	.	0
```

And now?

**Exercise**: What if we wanted to get a new bed file with only the neural crest specific enhancers? That is, we want to filter out anything that overlaps (by even 1 bp) with a H3K27Ac peak in stem cells. Look at the bedtools intersect -h output and see if you can figure out what option to use to do this.

Answer: `bedtools intersect -a data/human_enhancers.bed -b data/H1-hESC-H3K27Ac.bed -v > human_nconly_enhancers.bed`

Finally, we might want to generate a file where, for each neural crest enhancer in our list, we get the number of H1-hESC-H3K27Ac peaks it overlaps with. We'll use the `-c` option in bedtools for this. We'll also restrict overlaps to a reciprocal 20% -- so each feature in A has to overlap 25% of B, and vice versa, to count as an overlap.

```
bedtools intersect -a data/human_enhancers.bed -b data/H1-hESC-H3K27Ac.bed -f 0.20 -r -c > human_overlaps.txt
```

Now let's do the same thing with common enhancers:
```
bedtools intersect -a data/common_enhancers.bed -b data/H1-hESC-H3K27Ac.bed -f 0.20 -r -c > common_overlaps.txt
```

Get sequence from a bed file
--------------------

Another thing we might want to do is extract sequence data for particular intervals. You can do this in a couple of ways. By way of example we’ll use the dmel-genes.bed file. First, though, let’s pick 10 random genes to work with using bedtools sample, which selects random intervals from a file:

```
bedtools sample -i data/dmel-genes.bed -n 10 > dmel-10genes.bed
```

To get sequence for these genes, we can use bedtools getfasta:

```bedtools getfasta -fi data/dmel-all-chromosome-r6.20.fasta -bed dmel-10genes.bed -name > dmel-10genes-bedtools.fasta
```

We need to use the name option to include the bed name field in the fasta definition line.

We can also extract sequence with seqtk, which can take a bed file as an interval list:
```
seqtk subseq data/dmel-all-chromosome-r6.20.fasta dmel-10genes.bed > dmel-10genes-seqtk.fa
```

###maybe cut this section??
Anything different between the files?

We are going to look at one more use of bedtools, but first let’s do an exercise that will show you how you can string a lot of commands together to accomplish a task. The goal is to get the DNA sequence for 10 random exons from chromosome 2L in D. melanogaster. Take this in steps, and you can do testing using temporary files but in the end you’ll need to link everything in one command.

Shuffling Intervals
------------------

Let’s say we wanted to shuffle enhancers to random locations in the genome, but exclude them from known genes and introns (just allow them to fall in intergenic regions. Biologically this is not completely realistic, since enhancers can be in introns, but it will illustrate the procedure. We might want to do this to, for example, ask whether enhancers are closer to genes than might be expected by chance. To do this, we need to do a few steps:

First, we make a bed file of intergenic regions (note we need to include a sort command):

```
bedtools sort -i data/ucscGenes.bed -faidx data/hg19.genome | bedtools complement -i - -g data/hg19.genome > hg19_intergenic.bed
```

Next, we use shuffle with the -incl option to tell bedtools where it can put our shuffled elements:

```
bedtools shuffle -incl hg19_intergenic.bed -i data/common_enhancers.bed \
-noOverlapping -g data/hg19.genome > common_enhancers_random.bed
```
