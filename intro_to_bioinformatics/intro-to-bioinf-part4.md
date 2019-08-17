Read mapping
==============

We are going to cover how to map reads to a reference genome today. There are lots of reasons why you might want to do this (SNP calling, peak calling, RNA-seq expression estimates, etc). The pipeline we will go over is optimized for SNP calling, and will lead into tomorrow morning, when we'll discuss variant calling. But the basic ideas are similar for all data types.

We are going to be using a very small dataset downloaded from NCBI which represents resequencing data from 10 individual flycatchers, but filtered to include only a small region of the genome. Because of the way I generated this data, we are going to skip QC and trimming.

Set up
-------
We'll set up a directory structure to process each step of our pipline. Next, make a directory in your home dir or someone else you usually store files, and set it up as follows:

```
cd intro_bioinf_2019
mkdir -p read_mapping
cd read_mapping
mkdir -p 00_genome
mkdir -p 01_fastqs
mkdir -p 02_bams
mkdir -p 03_vcfs
```

This project organization is designed to keep things clean, and also remind you of the order of steps. Now we also need to copy over some data from elsewhere on Odyssey.

```
cp /n/holylfs/LABS/informatics/tsackton/workshops/read_mapping_variant_calling/genome/* 00_genome
cp /n/holylfs/LABS/informatics/tsackton/workshops/read_mapping_variant_calling/01_fastq/raw/* 01_fastqs
```

Mapping reads to a reference genome
----------------
Once you have your reads, you need to map them to a reference genome. There are many different aligners out there (e.g. BWA or bowtie2), but we recommend using BWA for SNP calling, so that read group information can be added during the alignment stage, without requiring a separate step with Picard’s AddOrReplaceReadGroups tool.

Before you can align your reads to a reference genome, you need to create an index. This only needs to be completed once per reference genome. BWA indexes are made from a FASTA genome file using bwa index:

```
cd 00_genome
bwa index -p ficAlb Falbicolis.chr5.fa
cd ..
```

The genome prefix should be a short identifier to be used as the prefix for all output files (e.g. prefix.bwt).

For most resequencing data, we want to use the bwa mem algorithm (for 70bp to 1Mbp query sequences) to map our reads. A typical command would be:

```
bwa mem -M -t 1 -R '@RG\tID:{FLOWCELL}.{LANE}\tPU:{FLOWCELL_BARCODE}.{LANE}.{SAMPLE}\tSM:{SAMPLE}\tPL:{PLATFORM}\tLB{LIBRARY}' <genome_prefix> <reads_1.fq> <reads_2.fq> > <samplename_bwa.sam>
```

There are many arguments available to use, as you can read in the manual. Some of the key arguments for these purposes are:

Argument	Description
-M	Mark shorter split hits as secondary - mandatory for Picard compatibility
-t <int>	Number of threads (default 1)
-R <str>	Read group information (see above for description)
-p	Specifies that fastq read 1 and read 2 files are interleaved, if only one fastq is specified and this command is not used, will assume single-end data

The output file format from BWA is a SAM (Sequence Alignment/Map) file format. This is a common file format, and detailed documentation can be found on the Samtools website. Samtools is part of a useful set of programs written to interact with high throughput sequencing data. The details of all you can do with this program are beyond the scope of this tutorial, but this program can be used to view, merge, calculate the depth of coverage, calculate other statistics, and index SAM-style files among other things.

For today, we’ll move to the 02_bams directory to do the mapping. We’ll start with a test run of a single sample, and then use a loop to do all of our input data at once. In real cases, you may use a loop to submit jobs with different parameters, or a job array, to do this efficiently.

```
cd 02_bams
bwa mem -M -t 1 -R '@RG\tID:COL_1\tSM:COL_1' ../00_genome/ficAlb ../01_fastqs/Falb_COL1.1.fastq ../01_fastqs/Falb_COL1.2.fastq \
> Falb_COL1.sam
```

Now, we are going to use a for loop to do this once for each sample. This is particularly efficient in this case since all the samples have numeric IDs, and writing a small loop in bash to iterate over a set of numbers is easy. The basic syntax is:

```
for INDEX in 1 2 3 4 5 31 32 33 34 35;
do
   echo $INDEX
done
```

When you run this, you should get a list of numbers printed to your screen. What is happening is that the INDEX variable is getting each value in the list in turn, and echo $VAR just prints what is in the variable $VAR.

Now, let’s do this for real:

```
for INDEX in 1 2 3 4 5 31 32 33 34 35;
do
   bwa mem -M -t 1 -R "@RG\tID:COL_$INDEX\tSM:COL_$INDEX" ../00_genome/ficAlb \
   ../01_fastqs/Falb_COL$INDEX.1.fastq \
   ../01_fastqs/Falb_COL$INDEX.2.fastq \
   > Falb_COL$INDEX.sam 2> Falb_COL$INDEX.log
done
```

Note that we are using the $INDEX variable in both the fastqs, the output, and the Read Group IDs.

You can examine some of the log files (e.g., with less) to see what is going on.

Sorting and Indexing
------------

Following alignment, you will need to sort the SAM file. We also recommend you store these types of alignment files in BAM format, which is similar to SAM format, but its binary equivalent, and therefore compressed and more efficient. As mentioned above, Samtools can be used to convert among file types, but we are working within the GATK pipeline for this tutorial, and so will work within the GATK/PicardTools universe.

You can automatically sort your SAM file by coordinate position (required for downstream analyses) and output the file in BAM format with PicardTools SortSam command.

```
picard SortSam \
I=samplename_bwa.sam \
O=samplename_sorted.bam \
SORT_ORDER=coordinate \
CREATE_INDEX=true
```

To use this BAM file, you also need to create a BAM index, so that software can efficiently access the compressed file. You notice that we automatically include index creation when sorting by specifying CREATE_INDEX=true. This is a universal command that can be applied to any PicardTools program. However, if you need to create an index separately, we do this with the BuildBamIndex command.

```
picard BuildBamIndex \
I=samplename_sorted.bam
```

Now let’s do this for our data. As before, we’ll run a test on one file, and then convert to a loop to process all BAMs.

```
picard SortSam \
I=Falb_COL1.sam \
O=Falb_COL1.sorted.bam \
SORT_ORDER=coordinate \
CREATE_INDEX=true
```

You can check to make sure you have a bam and a bam index (.bai) file with ls.

Now let’s do a loop:

```
for INDEX in 1 2 3 4 5 31 32 33 34 35
do
  picard SortSam \
  I=Falb_COL$INDEX.sam \
  O=Falb_COL$INDEX.sorted.bam \
  SORT_ORDER=coordinate \
  CREATE_INDEX=true
done
```

Alignment Metrics
-------

It may also be useful to calculate metrics on the aligned sequences. We can easily do this with the CollectAlignmentSummaryMetrics tool. Note that you can collect the alignment metrics on several different levels. In the below example, I’ve included metrics both at the sample and read group level. You also need to include the reference fasta file.

```
picard CollectAlignmentSummaryMetrics \
I=samplename_sorted.bam \
R=reference.fasta \
METRIC_ACCUMULATION_LEVEL=SAMPLE \
METRIC_ACCUMULATION_LEVEL=READ_GROUP \
O=samplename.alignment_metrics.txt
```

As an exercise, see if you can figure out how to modify the code block above to run on one of your samples. For an advanced challenge, use a loop to run on all your samples.

Deduplication
-----------

After alignment, sorting and indexing, it is necessary to identify any duplicate sequences from the same DNA fragment in your files that occur due to sample preparation (e.g. during PCR) or incorrect optical cluster identification during sequencing. This possibility is why it is important to identify read groups for different lanes of the same sample. This is also a useful point to merge together any BAM files from the same sample that are currently separated (demonstrated in example below). We identify duplicate sequences with MarkDuplicates, and additional details on how this is performed can be found in the tool documentation.

Note that it is not recommended to actually remove the duplicate sequences from the file, but simply to mark the flags appropriately in the BAM file, so that those sequences are ignored downstream. If using tools other than those we recommend here, make sure they can identify these flags. These sequences can also be removed later should the need arise.

At this point, we should be familiar with the loop command, so we are not going to run a test with just a single sample. Here is the full command:

```
for INDEX in 1 2 3 4 5 31 32 33 34 35
do
  picard MarkDuplicates \
  TMP_DIR=tmp \
  I=Falb_COL$INDEX.sorted.bam \
  O=Falb_COL$INDEX.dedup.bam \
  METRICS_FILE=Falb_COL$INDEX.dedup.metrics.txt \
  REMOVE_DUPLICATES=false \
  TAGGING_POLICY=All
done
```

We also recommend creating a deduplications metrics file, which will report the proportion and type of duplicate sequences in your sample and read groups.

Following deduplication make sure to sort and index your file again, as shown in the above section.

```
for INDEX in 1 2 3 4 5 31 32 33 34 35
do
  picard SortSam \
  I=Falb_COL$INDEX.dedup.bam \
  O=Falb_COL$INDEX.final.bam \
  SORT_ORDER=coordinate \
  CREATE_INDEX=true
done
```

Validating BAM files
--------

Once you are done with the above steps, it is best practice to validate your BAM file, to make sure there were not issues or mistakes associated with previous analyses. This is done with ValidateSamFile.

```
picard ValidateSamFile \
I=sample.dedup.sorted.bam \
O=sample.validate.txt \
MODE=SUMMARY
```

Again, in the interest of time we will skip this today.

At this point, we have an analysis-ready BAM. For some workflows, we would do indel realignment and base quality score recalibration at this point. These are described in more detail on the Informatics Pop Gen Tutorial, from which this workshop is derived.

With larger datasets, you would also want to remove some intermediate files once you have successfully validated your final bams.

Tomorrow, we will see how to use these bam files to call SNPs. For the rest of today, we are just going to look at a few additional things we can do with samtools to compute things about our BAM files.

Computing coverage
------
A common task is understanding the coverage of our sequencing reads across the genome.

```
samtools depth data/Falb_COL2.final.bam > Falb_COL2.depth
bedtools genomecov -ibam data/Falb_COL2.final.bam -bga -pc -g data/Falb_region.genome > Falb_COL2_depth.bedGraph
```

BAM to SAM conversion / BAM to BED conversion
----
One obvious task is to convert between sam and bam. We did this with Picard above, but we can also do it with samtools:

```
samtools view data/Falb_COL2.final.bam | head
samtools view -b data/Falb_COL2.final.sam > Falb_COL2.new.bam
```

We can also convert a bam file to a bed file, which can be useful for a number of reasons. We’ll do this in two steps: first, we need to produce a bam file that is sorted by query name (not start position), and excludes reads that are not properly paired:

```
samtools view -b -f 0x2 data/Falb_COL2.final.bam | samtools sort -n - | samtools fixmate -r - Falb_COL2.querysort.bam
```

Then we can use bamtobed to produce a bed file, where each interval is a read:

```
bedtools bamtobed -i Falb_COL2.querysort.bam -bedpe > Falb_COL2.bed
```
