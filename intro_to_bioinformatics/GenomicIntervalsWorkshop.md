Genomic ranges workshop
=========================

Many kinds of genomic data can be represented as **ranges**:
- Genes
- Chip-Seq peaks
- Promoters
- SNPs
- etc.

Essentially anything that represents an annotation on a genome sequence is probably best expressed as a range.
This means that working with range data is a very common problem in genomic analysis.

Today we are going to focus on a dataset recently published in [Cell](http://dx.doi.org/10.1016/j.cell.2015.08.036) that describes an extensive set of experiments to characterize human and chimpanzee neural crest cells and identify shared and species-specific enhancers. As part of their work, they provided a list of the genomic positions for 1000 human-specific enhancers, and 5000 enhancers shared between species, which will be the basis of our analysis.

Set up
------

These data are already on the cluster, so let's start by logging in to Odyssey and requesting an interactive node:


```
ssh tsackton@login.rc.fas.harvard.edu
srun -p interact --pty --mem 2000 -t 0-02:00 /bin/bash
```

Now let's set up the files:

```
cd /n/regal/informatics/workshops/GenomicIntervals_20151006/Users
mkdir tsackton
cp ../*.bed tsackton
cd tsackton
ls
```

You should have four files in your directory, all with the .bed extension. [BED files](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) are a standard format for storing range data in genomics. The minimal BED file is just a whitespace-delimited file with three columns (sequenceID, start position, and end position), where each row specifies an interval. The four files we have are three BED files with enhancer positions from various experiments, and one with UCSC genes. 

BED files can contain up to 9 additional optional fields, which have definitions specified by the order in which they appear (and if a higher number field appears, you have to include all the lower-numbered fields too). The 9 optional fields are:
- 4: name
- 5: score
- 6: strand
- 7: thickStart
- 8: thickEnd
- 9: itemRgb
- 10: blockCount
- 11: blockSizes
- 12: blockStarts

Fields 7, 8, and 9 reflect the origin of BED files in displaying information on a genome browser and are really only used for display purposes. Fields 10, 11, and 12 are used for specifying things like multi-exon genes in a single row of a BED file. 

> **Exercise**: Look at each bed file in your directory (use head or less) and figure out > which fields are included in each one.

**Note**: unlike GFF/GTF (another common file for storing ranges and annotation data in genomics), BED files represent ranges as zero-based, right open intervals (GFF files represent ranges as 1-based, closed intervals). So in a BED file the first 100 bases a chromosome would be represented by the interval 0-100 (which includes bases 0 through 99), whereas in a GFF file the same interval would be represented as 1-100 (which includes bases 1 through 100). *If you use specialized tools to manipulate these files, you don't have to worry about this complication, but if you write custom tools be careful about off-by-one errors.*

We will use two different specialized tools today: bedtools (a stand-alone program), and the GenomicRanges package in R/Bioconductor. We'll start with bedtools.

Bedtools
----------

First let's load bedtools on the cluster:

```
source new-modules.sh
module load bedtools2/2.25.0-fasrc01
```

Bedtools has a fairly comprehensive command line help, and also pretty good documentation [online](http://bedtools.readthedocs.org/en/latest/). Let's first just see what we can do by typing `bedtools` on the command line, which will give you a list of the sub-programs we can run. 

The first thing we'll do figure out what the closest gene is to each enhancer in our two enhancer bed files. **Exercise**: Scroll through the list of options and figure out what the best sub-command would be to do this.

To get help on a specific command, type `bedtools closest -h`

You'll see that for many commands, bedtools works by taking a 'query' bed file (specified with the -a option) and comparing it to a 'target' bed file (specified with the -b option). Although for some tools it doesn't matter which you put as -a and which as -b, in some cases it does matter so it is worth thinking through a bit.

**Exercise**: to find which gene is closest to each enhancer in a given BED file, what should be specify as -a and what as -b?

Answer:
```
bedtools closest -a human_enhancers.bed -b ucscGenes.bed | head

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

Answer: `bedtools closest -a human_enhancers.bed -b ucscGenes.bed | awk -v OFS='\t' '{print $1, $2, $3, $9, $5}' | head`

Okay, now let's make a file that contains the closest gene for each of the human-specific enhancers, in proper bed format, but we'll add two additional options (refer back to bedtools closest -h for the full option list): -d to give us the distance to the nearest gene, which we'll put in the 'score' field of the output bed file, and -t "first" means that if there are ties, we'll just keep one (so that we have each enhancer assigned to a single gene). 

`bedtools closest -a human_enhancers.bed -b ucscGenes.bed -d -t "first" | awk -v OFS='\t' '{print $1, $2, $3, $9, $10}' > hse_closest_genes.bed`

**Exercise**: What do we need to change to get the equivalent set but for the common enhancers?

Answer: The input file name for bed file 'A', the output file name, and the field numbers for the gene and score.
```
bedtools closest -a common_enhancers.bed -b ucscGenes.bed -d -t "first" | awk -v OFS='\t' '{print $1, $2, $3, $7, $8}' > ce_closest_genes.bed
```

Another common task is finding overlaps between two sets of intervals. For example, we might want to know which of our enhancers identified in neural crest cells are shared with enhancers identified in some other cell type. In your directory you should see another bed file, H1-hESC-H3K27Ac.bed, which has H3K27Ac peaks from human embryonic stem cells. We'll use bedtools to get the overlaps between this file and each of our enhancer files. 

The tool we'll use is intersect. Let's start by looking at the options: ```bedtools intersect -h```

The output of this can be a bit complicated, so let's start by just looking at a few things:

```
bedtools intersect -a human_enhancers.bed -b H1-hESC-H3K27Ac.bed | head -n 3

chr1	1590473	1590673	HSE853	17.3206898329
chr1	6336418	6336624	HSE394	14.8892086906
chr1	7404594	7404794	HSE315	24.228051023
```

What are we seeing here?

```
bedtools intersect -a human_enhancers.bed -b H1-hESC-H3K27Ac.bed -wa | head -n 3

chr1	1590473	1590673	HSE853	17.3206898329
chr1	6336418	6336624	HSE394	14.8892086906
chr1	7404594	7404794	HSE315	24.228051023

bedtools intersect -a human_enhancers.bed -b H1-hESC-H3K27Ac.bed -wb | head -n 3
chr1	1590473	1590673	HSE853	17.3206898329	chr1	1590264	1590895	.	1000	.
chr1	6336418	6336624	HSE394	14.8892086906	chr1	6160044	6490440	.	291	.
chr1	7404594	7404794	HSE315	24.228051023	chr1	7400919	7412366	.	279	.
```

What about now?

```
bedtools intersect -a human_enhancers.bed -b H1-hESC-H3K27Ac.bed -loj | head -n 3
chr1	1015066	1015266	HSE897	9.5706324138	.	-1	-1	.	-1	.
chr1	1590473	1590673	HSE853	17.3206898329	chr1	1590264	1590895	.	1000	.
chr1	2120861	2121064	HSE86	66.0424471614	.	-1	-1	.	-1	.

bedtools intersect -a human_enhancers.bed -b H1-hESC-H3K27Ac.bed -wao | head -n 3
chr1	1015066	1015266	HSE897	9.5706324138	.	-1	-1	.	-1	.	0
chr1	1590473	1590673	HSE853	17.3206898329	chr1	1590264	1590895	.	1000	.	200
chr1	2120861	2121064	HSE86	66.0424471614	.	-1	-1	.	-1	.	0
```

And now?

**Exercise**: What if we wanted to get a new bed file with only the neural crest specific enhancers? That is, we want to filter out anything that overlaps (by even 1 bp) with a H3K27Ac peak in stem cells. Look at the bedtools intersect -h output and see if you can figure out what option to use to do this.

Answer: `bedtools intersect -a human_enhancers.bed -b H1-hESC-H3K27Ac.bed -v > human_nconly_enhancers.bed`

Finally, we might want to generate a file where, for each neural crest enhancer in our list, we get the number of H1-hESC-H3K27Ac peaks it overlaps with. We'll use the `-c` option in bedtools for this. We'll also restrict overlaps to a reciprocal 20% -- so each feature in A has to overlap 25% of B, and vice versa, to count as an overlap. 

```
bedtools intersect -a human_enhancers.bed -b H1-hESC-H3K27Ac.bed -f 0.20 -r -c > human_overlaps.txt
```

Now let's do the same thing with common enhancers:
```
bedtools intersect -a common_enhancers.bed -b H1-hESC-H3K27Ac.bed -f 0.20 -r -c > common_overlaps.txt
```

There is a lot more that one can do with bedtools, but hopefully this serves as a useful introduction. Any specific questions before we move on to R?

Genomic Ranges in R/Bioconductor
--------------

Working with genomic ranges in bedtools is very powerful, especially as one can take advantage of UNIX pipes to link together analysis chains. But there are also some advantages to R, so we're going to do essentially the same analysis in R, on your laptop. First, close your interactive session, and then make a working directory somewhere on your computer, cd to it, and download all the files in the directory you just created to that working directory. 

```
cd /Users/tim/Projects/RC/workshops/20151006-GenomicIntervals
scp tsackton@login.rc.fas.harvard.edu:/n/regal/informatics/workshops/GenomicIntervals_20151006/Users/tsackton/*.* .
```

R has a whole set of packages specifically to deal with genomic ranges. We'll just scratch the surface today. Some other good introductions include:
- [Vince Buffalo's Genomic Ranges introduction](https://github.com/vsbuffalo/genomicranges-intro/blob/master/notes.Rmd)
- [The Genomic Ranges HOWTO](http://bioconductor.org/packages/release/bioc/vignettes/GenomicRanges/inst/doc/GenomicRangesHOWTOs.pdf)

Start by opening RStudio, and changing your working directory to the place where you put the files downloaded from the cluster. Now we need to install the necessary R packages from Bioconductor:

```
setwd("~/Projects/RC/workshops/20151006-GenomicIntervals")
source("http://bioconductor.org/biocLite.R")
biocLite("GenomicRanges")
biocLite("rtracklayer")
biocLite("TxDb.Hsapiens.UCSC.hg19.knownGene")
```

And load them (loading rtracklayer will load genomic ranges and related packages as dependencies:
```
library(rtracklayer)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
```

rtracklayer provides an interface for loading BED files into R, into a special data class called a GRanges object. GRanges objects have methods that let you access positional information directly, and also allow operations on ranges such as overlaps. They are extremely powerful and we are just barely scratching the surface today. They are built upon a more general class called 'IRanges'. 

```
humEn<-import("human_enhancers.bed") 
sharedEn<-import("common_enhancers.bed")
```

Let's start by exploring these objects a bit. Just typing the name of an object will give us the first few and last few rows:
```
humEn

GRanges object with 1000 ranges and 2 metadata columns:
         seqnames                 ranges strand   |        name     score
            <Rle>              <IRanges>  <Rle>   | <character> <numeric>
     [1]     chr1     [1015067, 1015266]      *   |      HSE897  9.570632
     [2]     chr1     [1590474, 1590673]      *   |      HSE853 17.320690
     [3]     chr1     [2120862, 2121064]      *   |       HSE86 66.042447
     [4]     chr1     [6336419, 6336624]      *   |      HSE394 14.889209
     [5]     chr1     [7404595, 7404794]      *   |      HSE315 24.228051
     ...      ...                    ...    ... ...         ...       ...
   [996]     chrX [114733852, 114734051]      *   |      HSE986  7.992718
   [997]     chrX [124903683, 124903882]      *   |      HSE436 12.664842
   [998]     chrX [128211487, 128211686]      *   |      HSE815 14.361143
   [999]     chrX [128270782, 128270981]      *   |      HSE952 10.135859
  [1000]     chrX [130964621, 130964822]      *   |      HSE980 11.910385
  -------
  seqinfo: 23 sequences from an unspecified genome; no seqlengths
```

We can also access information with accessor functions:
```
head(score(humEn))
[1]  9.570632 17.320690 66.042447 14.889209 24.228051 17.519263
head(width(humEn))
[1] 200 200 203 206 200 200
head(start(humEn))
[1]  1015067  1590474  2120862  6336419  7404595 11941326
```

This is useful for making plots and other summaries:
```
plot(density(score(humEn)))
hist(width(humEn))
table(seqnames(humEn))
```

**Exercise**: What proportion of human-specific enhancers are on the X chromosome? What about shared enhancers? Are these proportions statistically different?

Answer:
```
prop.table(table(seqnames(sharedEn) == "chrX"))

FALSE  TRUE 
0.979 0.021 

prop.table(table(seqnames(humEn) == "chrX"))

FALSE  TRUE 
0.979 0.021 

fisher.test(matrix(c(table(seqnames(sharedEn) == "chrX"),table(seqnames(humEn) == "chrX")),nrow=2,ncol=2))
	Fisher's Exact Test for Count Data

data:  
p-value = 1
alternative hypothesis: true odds ratio is not equal to 1
95 percent confidence interval:
 0.5912779 1.6182666
sample estimates:
odds ratio 
         1 
```

We can also subset genomic ranges like normal R objects. For example, if we want the human-specific neural crest enhancers on chromsome 1:
```
humEn[seqnames(humEn)=="chr1",]

GRanges object with 83 ranges and 2 metadata columns:
       seqnames                 ranges strand   |        name     score
          <Rle>              <IRanges>  <Rle>   | <character> <numeric>
   [1]     chr1     [1015067, 1015266]      *   |      HSE897  9.570632
   [2]     chr1     [1590474, 1590673]      *   |      HSE853 17.320690
   [3]     chr1     [2120862, 2121064]      *   |       HSE86 66.042447
   [4]     chr1     [6336419, 6336624]      *   |      HSE394 14.889209
   [5]     chr1     [7404595, 7404794]      *   |      HSE315 24.228051
   ...      ...                    ...    ... ...         ...       ...
  [79]     chr1 [242513808, 242514007]      *   |       HSE29 41.864605
  [80]     chr1 [243428177, 243428376]      *   |      HSE909  8.976846
  [81]     chr1 [245362207, 245362406]      *   |      HSE514 20.851024
  [82]     chr1 [245748623, 245748822]      *   |      HSE400 15.880901
  [83]     chr1 [245751246, 245751445]      *   |      HSE614 11.997767
  -------
  seqinfo: 23 sequences from an unspecified genome; no seqlengths
```


Now, let's replicate our bedtools analysis:

To get the closest gene to each enhancer, we first need to get the location of genes. We could do this by reading in a bed file, but within R there are easier ways. In particular, for humans and a few other commonly used organisms, Bioconductor maintains genomic annotation packages, and we can use accessor functions to get genomic ranges objects from them. This is a whole separate workshop, so I'm just going to give you the code:

```
hg<-transcripts(TxDb.Hsapiens.UCSC.hg19.knownGene)
```

**Exercise**: Make a histogram of the lengths of human transcripts.

Answer:
```
hist(width(hg))
```

**Exercise**: What is the ID of the longest transcript in humans?

Answer:
```
hg[width(hg)==max(width(hg)),]
```

Okay, so the goal here is actually to find the closest gene to each human-specific enhancer. We use the nearest function to do this:
```
humEn.trans<-nearest(humEn, hg)
```

This will give us just a list of the row ids in hg that are nearest to each row of the humEn object, but what if we want something more useful?

```
humEn.trans.df<-data.frame(en=humEn$name, trans=hg[humEn.trans,]$tx_name)
humEn.trans.gr<-hg[humEn.trans,]
humEn$tx_name=humEn.trans.gr$tx_name
```

Finally, let's get the intersection between our neural crest enhancers and the embryonic stem cell H3K27Ac peaks. First we need to read in the hESC bed file:
```
hESC<-import("H1-hESC-H3K27Ac.bed")
```

Finding overlaps uses, perhaps unsurprisingly, the findOverlaps function:
```
hESC.overlap<-findOverlaps(humEn,hESC)
```

To get counts we can use the countOverlaps function:
```
hESC.overlap.counts<-countOverlaps(humEn,hESC)
hESC.shared.counts<-countOverlaps(sharedEn,hESC)
```

Are human-specific neural crest enhancers more or less likely to overlap stem cell peaks than common neural crest enhancers?
```
fisher.test(cbind(table(hESC.overlap.counts>=1),table(hESC.shared.counts>=1)))
```

**Exercise**: Read in your bedtools count output and compare to the R output. Are there any differences? If so why do you think that is?

