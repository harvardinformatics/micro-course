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
```

**Exercise**: in this case, does bedtools output a bed file? How could we convert this to a proper bed file, where we keep the score but replace the names (HSE####) with the closest gene? (Remember `awk '{print $1}` prints the first field of a file)?

Answer:
```
```

Okay, now let's make a file that contains the closest gene for each of the human-specific enhancers, in proper bed format, but we'll add two additional options (refer back to bedtools closest -h for the full option list): -d to give us the distance to the nearest gene, which we'll put in the 'score' field of the output bed file, and -t "first" means that if there are ties, we'll just keep one (so that we have each enhancer assigned to a single gene).

`bedtools closest -a data/human_enhancers.bed -b data/ucscGenes.bed -d -t "first" | awk -v OFS='\t' '{print $1, $2, $3, $9, $10}' > hse_closest_genes.bed`

**Exercise**: What do we need to change to get the equivalent set but for the common enhancers? Run bedtools first to see what the output looks like (piping to less or head), and then edit our awk command to get the correct fields to produce a bed file.

Answer:
```
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

Answer:
```
```


Finally, we might want to generate a file where, for each neural crest enhancer in our list, we get the number of H1-hESC-H3K27Ac peaks it overlaps with. We'll use the `-c` option in bedtools for this. We'll also restrict overlaps to a reciprocal 20% -- so each feature in A has to overlap 25% of B, and vice versa, to count as an overlap.

```
bedtools intersect -a data/human_enhancers.bed -b data/H1-hESC-H3K27Ac.bed -f 0.20 -r -c > human_overlaps.txt
```

Now let's do the same thing with common enhancers:

```
```

Get sequence from a bed file
--------------------

Another thing we might want to do is extract sequence data for particular intervals. You can do this in a couple of ways. By way of example we’ll use the dmel-genes.bed file. First, though, let’s pick 10 random genes to work with using bedtools sample, which selects random intervals from a file:

```
bedtools sample -i data/dmel-genes.bed -n 10 > dmel-10genes.bed
```

To get sequence for these genes, we can use bedtools getfasta:

```
bedtools getfasta -fi data/dmel-all-chromosome-r6.20.fasta -bed dmel-10genes.bed -name > dmel-10genes-bedtools.fasta
```

We need to use the name option to include the bed name field in the fasta definition line.

We can also extract sequence with seqtk, which can take a bed file as an interval list:

```
seqtk subseq data/dmel-all-chromosome-r6.20.fasta dmel-10genes.bed > dmel-10genes-seqtk.fa
```

Anything different between the files?

We are going to look at one more use of bedtools, but first let’s do an exercise that will show you how you can string a lot of commands together to accomplish a task. The goal is to get the DNA sequence for 10 random exons from chromosome 2L in D. melanogaster. Take this in steps, and you can do testing using temporary files but in the end you’ll need to link everything in one command.

Shuffling Intervals
------------------

Let’s say we wanted to shuffle enhancers to random locations in the genome, but exclude them from known genes and introns (just allow them to fall in intergenic regions. Biologically this is not completely realistic, since enhancers can be in introns, but it will illustrate the procedure. We might want to do this to, for example, ask whether enhancers are closer to genes than might be expected by chance. To do this, we need to do a few steps:

```
bedtools complement -i data/ucscGenes.bed -g data/hg19.genome
```

There is an error; we need to sort the bed file. We can use bedtools sort and pipe the output to the input of bedtools complement.
Note that when we pipe to bedtools, we need to use `-` to indicate which file is coming via the previous command.
E.g., ```bedtools complement -i -``` would indicate that the input to complement is taken from stdout of the piped command.

Now combining:
```
bedtools sort -i data/ucscGenes.bed -faidx data/hg19.genome | bedtools complement -i - -g data/hg19.genome > hg19_intergenic.bed
```

Next, we use shuffle with the -incl option to tell bedtools where it can put our shuffled elements:

```
bedtools shuffle -incl hg19_intergenic.bed -i data/common_enhancers.bed \
-noOverlapping -g data/hg19.genome > common_enhancers_random.bed
```

We'll end there. After lunch we'll pick up with read mapping and a little bit of information on working with bam and sam files.
