## Table of Contents
[Whole-genome Pop Gen Sequencing Overview](#overview)<br>
[Experimental Design](#design)<br>
[Sequence Reads](#reads)<br>
[Quality Control](#qc)<br>
[Preprocessing](#preprocess)<br>
[Indel Realignment](#realignment)<br>
[Base Quality Score Recalibration](#bqsr)<br>
[Variant Calling](#variantcalling)<br>
[Data Filtering](#filtering)<br>
[Next Steps](#next)<br>
[References](#references)<br>

## Whole-genome resequencing population genomics overview <a name="overview"></a>

Population genetics can be used to identify genetic variation within and between populations, and with DNA sequencing becoming less expensive, more researchers are turning to whole-genome resequencing to understand genome-wide variation. The objective of this tutorial is to familiarize users with the process of obtaining analysis-ready VCF files from population genomic whole-genome resequencing data. The tutorial is based on the [GATK's best practices pipeline for Germline SNP and Indel Discovery](https://software.broadinstitute.org/gatk/best-practices/bp_3step.php?case=GermShortWGS), however, geared toward non-human organisms. We also address low-coverage whole-genome resequencing data in the tutorial, as we expect this data type to be common for our users. In addition to [fastq sequencing data files](#reads), it is also necessary to have a reference genome fasta file for this pipeline. If the reference exists but you don't have it in hand, you can download the fasta file from [that organism's genome page from NCBI](https://www.ncbi.nlm.nih.gov/genome/browse/#!/overview/).

## Experimental design <a name="design"></a>

There are numerous preparation methods, (e.g. [Nextera](https://www.illumina.com/products/by-type/sequencing-kits/library-prep-kits/nextera-dna.html), [Kapa](http://sequencing.roche.com/en/products-solutions/by-category/library-preparation/dna-library-preparation.html) ) to construct sequencing libraries from DNA extraction. These laboratory methods are beyond the scope of this tutorial. However, we will address a few aspects of study design when designing an experiment.

### 1. Pooled sequencing vs. individually barcoding samples

One decision researchers need to make when designing their resequencing experiments is whether to pool unbarcoded individuals together in a single sequencing library (termed **Pool-seq**), or to individually barcode each individual, enabling researchers to demultiplex these individuals for downstream analyses even if they are pooled for sequencing itself. There are pros and cons to both approaches, but essentially the decision comes down to cost and research objective. Many papers have been written about the pros and cons of pooled sequencing, and [Schlötterer et al. 2014](https://www.nature.com/articles/nrg3803) provides a nice review and comparison with other methods. Briefly, we outline some of the pros and cons below:

* **Pooled sequencing**: The main advantage for this approach is cost savings on library preparation. If large sample sizes are required for the research objectives, library preparation costs can quickly become a limiting factor. By pooling large numbers of individuals in a single population, researchers would only need to prepare a single sequencing library per pool. However, this method has limitations on possible downstream analyses and potential sequencing biases. This method can yield estimates of allele frequencies from a pooled population, but few statistics beyond that (e.g. haplotype information, linkage disequilibrium). Pool-seq also works best when large numbers of individuals (>40) are pooled together, with pools on the order of hundreds or thousands of individuals being ideal. One of the biggest drawbacks of Pool-seq is that unequal individual representation will produce biases in allele frequency estimates, and without barcodes it is impossible to know if this has occurred. This is less likely to happen with larger sample sizes.

* **Individually barcoded sequencing**: The main advantage of this approach is that individually barcoding reads means that variants can be called for individuals, and with sufficient coverage (see below), it is possible to obtain haplotype information or other useful statistics. As mentioned above, the main drawback of this method is the cost of library preparation. This is increasingly less expensive however, either because of new kits becoming available, or the ability to split library preparation reagents into multiple microreactions (e.g. [Baym et al. 2015](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4441430/)). Therefore, in the tutorial, we focus on methods for creating VCFs from individually barcoded samples.

### 2. Sample sizes

Determining how many individuals you need to sequence depends on what types of analysis you wish to conduct downstream. If the objective of the study is to describe population structure and genetic diversity, very few individuals per population are needed because whole-genome sequencing provides so much information per individual (e.g. [Nazareno et al. 2017](http://onlinelibrary.wiley.com/doi/10.1111/1755-0998.12654/full)). If the goal of the study is to perform detailed demographic inference (e.g. with the site frequency spectrum via dadi or fastsimcoal2), small numbers of individuals may be sufficient for detecting older events or testing different models, but larger numbers of individuals may be necessary to detect recent events or estimate parameters (e.g. [Robinson et al. 2014](https://bmcevolbiol.biomedcentral.com/articles/10.1186/s12862-014-0254-4)).

Identifying allele frequency shifts at specific sites (e.g. looking for Fst outliers) or GWAS-types of analyses require larger population sizes to have enough power to detect significant differences with FDR corrections for millions of sites.

### 3. Sequencing depth

Ideally, to confidently call variants from whole genome resequencing data, diploid organisms should be sequenced to 30x coverage. However, due to limited budgets and different study goals, it is often possible to sequence to much lower coverage. For many population genomic goals, given a set amount of sequencing to work with (e.g. 100x coverage of the target genome), it is often more advantageous to sequence more individuals (100 individuals to 1x coverage), to more accurately infer population genetic parameters ([Buerkle and Gompert 2012](http://onlinelibrary.wiley.com/doi/10.1111/mec.12105/abstract)). Because of limited budgets despite falling sequencing costs, an increasing number of tools are available to make use of low-coverage whole-genome resequencing for population genomic inference. For example, the [ANGSD](http://www.popgen.dk/angsd/index.php/ANGSD) and [NGSTools](https://github.com/mfumagalli/ngsTools) packages allow one to calculate site frequency spectra, diversity statistics, PCA, and admixture analysis among others based entirely on genotype likelihood scores. Other packages, like [MAPGD](https://github.com/LynchLab/MAPGD), allow one to calculate linkage disequilibrium and relatedness using genotype likelihoods. By not actually calling genotypes, and instead inferring parameters from genotype likelihoods across individuals in a population, these programs avoid many of the biases associated with low-coverage genome data (e.g. [Han et al. 2014](https://academic.oup.com/mbe/article/31/3/723/1007998)).

### 4. Sequencing mode

For whole-genome resequencing studies, it is almost always recommended to use paired-end sequencing. As genome coverage is generally a limiting factor, the cost per base is much less for paired-end than single-end data. In addition, paired-end data generally provide better abilities to map reads to reference genomes, which is highly advantageous, especially for low-coverage data.

## Quality control <a name="qc"></a>

Before beginning analysis, it is a good idea to assess the general quality of the raw sequence data. One fast and easy to use program for this is [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/). If you sequence your data at the Bauer core, this will be automatically run for you. 
    
This will produce an HTML file as output, which can be viewed in any web browser. It will provide a number of useful summary statistics, including graphical representations of your data. One of the most useful metrics in this context is to see how the quality of your sequence data varies along the length of the read. It is normal to see quality drop along a sequence, but large drops in the "red" zone, especially for older sequence data, may indicate that uniform trimming to exclude those low-quality bases is a good idea. See more about trimming below.

## Running the GATK/PicardTools Pipeline on Odyssey

##### The commands and information given below are specific to GATK and PicardTools version 3. These programs have been rolled together in GATK4, and while the commands will be similar, there will be some changes in the command structure. Although we have started using GATK4 in production work, many workflows continue to use GATK3 and for teaching purposes we are continuing with GATK3. If you are interested in updating your workflow to GATK4 and are stuck, we can help. If you download the latest version of GATK, you will get GATK4. You can download previous versions of GATK (e.g. GATK 3.8 [here](https://software.broadinstitute.org/gatk/download/archive)).

A few notes on running **GATK** and **PicardTools** commands on Odyssey. 

**GATK** and **PicardTools** are built with java, and so when running the *jar* file (e.g. `java -jar picard.jar <PicardTool>`), you can include a few extra [options](https://docs.oracle.com/javase/7/docs/technotes/tools/windows/java.html) to pass to java that are especially applicable to running these programs on Odyssey. First, you can add a memory limit to java, for example requiring java to use no more than 4GB memory: `-Xmx4g`. This can help ensure your program does not use more memory than you request, resulting in job termination on Odyssey. 

The second is: `-XX:ParallelGCThreads=1`. This command limits the number of java "Garbage Collector" threads running for each task. Based on what you set this number to be (we recommend 1 or 2), you should make sure to request that that number +1 for the number of cores you request in your submission script (e.g. `-n 2`). If you do not request this, java will start using many threads, and may cause your script to unexpectedly fail.

Finally, you must be using java version 8 (or development kit 1.8). So, be sure to load a java module before running any commands. While there are GATK modules installed on Odyssey, it is simple to download the latest versions yourself. Just be sure to download and specify the picard.jar and GenomeAnalysisTK.jar files when running commands.

An example PicardTools command with these two variables is: 

    module load jdk/1.8.0_45-fasrc01
    java -Xmx4g -XX:ParallelGCThreads=1 -jar ~/path/to/picard.jar <PicardTool> \
    I=<Input> \
    1=<Option1> \
    2=<Option2> \
    O=<Output>
    _
Note that the options associated with the program can all go in a single line without the `\` at the end, but it makes a script easier to read and interpret to break the command across lines. Commands that require use input are demonstrated by `<>`. You should fill in these options with your own filenames or commands.

## Preprocessing <a name="preprocess"></a>

Preprocessing sequencing reads is much simpler when dealing with resequencing data compared to building a genome. Most resequencing libraries will be created with large enough insert sizes so that paired-end reads will not overlap, and adapter contamination is not an issue (however, if this is not the case with your data, please visit the options for adapter removal in the [ATAC-seq workflow](https://informatics.fas.harvard.edu/peak-calling-workflow.html#qc)). 

Note that if for these reasons you also have many overlapping reads, in addition to trimming adapters, you also might want to merge the reads into single fragments. If you plan to call variants with a tool like **HaplotypeCaller** from the GATK pipeline as we describe below, this is not necessary as the program will account for overlapping bases and not inflate the variant quality scores. However, if the program you plan to use does not account for overlapping bases (e.g. **ANGSD**) you may want to use software like **NGmerge** in *stitch mode*, created by the Informatics group, to merge the reads. 

    module load NGmerge
    NGmerge -1 <fastq_read_1> -2 <fastq_read_2> -o <output_fasta_name>

See the documentation for NGmerge `NGmerge -h` for additional parameter options, such as the number of threads to use.

Note that additional trimming with resequencing data is not usually necessary, as many variant callers (e.g. **HaplotypeCaller**) take quality scores into account. Others (e.g. **ANGSD**), can trim reads during the data filtering step. For that reason, we do not recommend trimming here.

### Read Group Identifiers

Read group identifiers are used to identify sequence data by sequencing technology (e.g. Illumina), flow cell, lane, sample ID, and library. Using these identifiers ensures that batch effects, or biases in the data that might have been introduced at different stages of the sequencing process can be properly accounted for. Detailed documentation on Read Groups can be found on the [GATK website](https://gatkforums.broadinstitute.org/gatk/discussion/6472/read-groups).

The most common and recommended read groups are:

* `ID` : **Read Group Identifier**<br>
    A unique identifier for each read group. The convention for Illumina data is {FLOWCELL}.{LANE}.

* `PU`: **Platform Unit**<br>
    A sample/library specific identifier, specified with: {FLOWCELL_BARCODE}.{LANE}.{SAMPLE}. The *flowcell barcode* is a unique identifier for a flow cell, *lane* is the lane of that flowcell, and *sample* is the sample or library specific identifier.

* `SM`: **Sample**<br>
    The name of the sample represented by this read group. This will be the name used in the sample column of the VCF file.

* `PL`: **Platform**<br>
    The sequencing technology used to create the data. Current valid values: ILLUMINA, SOLID, LS454, HELICOS, and PACBIO.

* `LB`: **Data Preparation Library Identifier**<br>
    The library preparation identifier. This is used by MarkDuplicates to identify which read groups contain molecular (e.g. PCR) duplicates.

The read group information can be found in the file header (look for `@RG`) and the `RG:Z` tag for each sequence record. This information is not automatically added to Fastq files following sequencing, but needs to be added either when mapping with **BWA** or separately after mapping with Picard's [AddOrReplaceReadGroups](http://broadinstitute.github.io/picard/command-line-overview.html#AddOrReplaceReadGroups) tool.
 
If you don't know the information on the flowcell and lane for your data, you can derive the information from the sequence headers found in a Fastq file, as described in the [Sequence reads](#reads) section.

##Pipeline for today's workshop

For today, we will start here at read mapping. We are going to be using a very small dataset downloaded from NCBI which represents resequencing data from 10 individual flycatchers, but filtered to include only a small region of the genome. Because of the way I generated this data, we are going to skip QC and trimming.

###Set up

To get set up, first we will need to log into Odyssey and request an interactive node. For larger scale projects, we'd be doing this with sbatch scripts but for today most steps will take less than a minute to run. 

    srun -p test --pty --mem 4000 -n 2 -N 1 -t 0-02:00 /bin/bash

This will request an interactive node with 4GB of memory and 2 cores, for 2 hours. 

Next, make a directory in your home dir or someone else you usually store files, and set it up as follows:

    cd workshops
    mkdir -p rmvc_fall2018
    cd rmvc_fall2018
    mkdir -p 00_genome
    mkdir -p 01_fastqs
    mkdir -p 02_bams
    mkdir -p 03_vcfs

This project organization is designed to keep things clean, and also remind you of the order of steps. Now we also need to copy over some data from elsewhere on Odyssey.

    cp /n/holylfs/LABS/informatics/workshops/read_mapping_variant_calling/genome/* 00_genome
    cp /n/holylfs/LABS/informatics/workshops/read_mapping_variant_calling/01_fastq/raw/* 01_fastqs

##Trimming

module load NGmerge/0.2-fasrc01
NGmerge -1 Falb_COL1.1.fastq -2 Falb_COL1.2.fastq -a -v -o Falb_COL1.trimmed

for INDEX in 1 2 3 4 5 31 32 33 34 35;
do
	NGmerge -1 Falb_COL$INDEX.1.fastq -2 Falb_COL$INDEX.2.fastq -a -v -o Falb_COL$INDEX.trimmed 2> COL_$INDEX.log
done
N

### Mapping reads to a reference genome

Once you have your reads, you need to map them to a reference genome. There are many different aligners out there (e.g. [BWA](http://bio-bwa.sourceforge.net) or [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)), but we recommend using **BWA** so that read group information can be added during the alignment stage, without requiring a separate step with Picard's [AddOrReplaceReadGroups](http://broadinstitute.github.io/picard/command-line-overview.html#AddOrReplaceReadGroups) tool.

Before you can align your reads to a reference genome, you need to create an index. This only needs to be completed once per reference genome. BWA indexes are made from a FASTA genome file using `bwa index`:

    module load bwa/0.7.15-fasrc02
    cd 00_genome
    bwa index -p ficAlb Falbicolis.chr5.fa
    cd ..

The genome prefix should be a short identifier to be used as the prefix for all output files (e.g. `prefix.bwt`).

For most resequencing data, we want to use the `bwa mem` algorithm (for 70bp to 1Mbp query sequences) to map our reads. A typical command would be:

    bwa mem -M -t 1 -R '@RG\tID:{FLOWCELL}.{LANE}\tPU:{FLOWCELL_BARCODE}.{LANE}.{SAMPLE}\tSM:{SAMPLE}\tPL:{PLATFORM}\tLB{LIBRARY}' <genome_prefix> <reads_1.fq> <reads_2.fq> > <samplename_bwa.sam>

There are many arguments available to use, as you can read in the [manual](http://bio-bwa.sourceforge.net/bwa.shtml). Some of the key arguments for these purposes are:

| Argument   | Description                                  |
|:----------:|----------------------------------------------|
| `-M`       | Mark shorter split hits as secondary - mandatory for Picard compatibility |
| `-t <int>` | Number of threads (default 1) |
| `-R <str>` | Read group information (see above for description) |
| `-p`       | Specifies that fastq read 1 and read 2 files are interleaved, if only one fastq is specified and this command is not used, will assume single-end data |

The output file format from BWA is a SAM (Sequence Alignment/Map) file format. This is a common file format, and [detailed documentation](https://samtools.github.io/hts-specs/SAMv1.pdf) can be found on the Samtools website. **[Samtools](http://www.htslib.org)** is part of a useful set of programs written to interact with high throughput sequencing data. The details of all you can do with this program are beyond the scope of this tutorial, but this program can be used to view, merge, calculate the depth of coverage, calculate other statistics, and index SAM-style files among other things.

For today, we'll move to the 02_bams directory to do the mapping. We'll start with a test run of a single sample, and then use a loop to do all of our input data at once. In real cases, you may use a loop to submit jobs with different parameters, or a job array, to do this efficiently.

    cd 02_bams
    bwa mem -M -t 1 -R '@RG\tID:COL_1\tSM:COL_1' ../00_genome/ficAlb ../01_fastqs/Falb_COL1.trimmed_1.fastq ../01_fastqs/Falb_COL1.trimmed_2.fastq \
    > Falb_COL1.sam
    
Now, we are going to use a for loop to do this once for each sample. This is particularly efficient in this case since all the samples have numeric IDs, and writing a small loop in bash to iterate over a set of numbers is easy. The basic syntax is:

    for INDEX in 1 2 3 4 5 31 32 33 34 35;
    do
       echo $INDEX
    done

When you run this, you should get a list of numbers printed to your screen. What is happening is that the INDEX variable is getting each value in the list in turn, and echo $VAR just prints what is in the variable $VAR.

Now, let's do this for real:

    for INDEX in 1 2 3 4 5 31 32 33 34 35;
    do
       bwa mem -M -t 1 -R "@RG\tID:COL_$INDEX\tSM:COL_$INDEX" ../00_genome/ficAlb \
       ../01_fastqs/Falb_COL$INDEX.trimmed_1.fastq \
       ../01_fastqs/Falb_COL$INDEX.trimmed_2.fastq \
       > Falb_COL$INDEX.sam 2> Falb_COL$INDEX.log
    done

Note that we are using the $INDEX variable in both the fastqs, the output, and the Read Group IDs.

You can examine some of the log files (e.g., with less) to see what is going on.

### Sorting and Indexing

Following alignment, you will need to sort the SAM file. We also recommend you store these types of alignment files in BAM format, which is similar to SAM format, but its binary equivalent, and therefore compressed and more efficient. As mentioned above, **[Samtools](http://www.htslib.org)** can be used to convert among file types, but we are working within the GATK pipeline for this tutorial, and so will work within the GATK/PicardTools universe. 

You can automatically sort your SAM file by coordinate position (required for downstream analyses) and output the file in BAM format with PicardTools **[SortSam](http://broadinstitute.github.io/picard/command-line-overview.html#SortSam)** command. 

    java -jar ~/path/to/picard.jar SortSam \
    I=samplename_bwa.sam \
    O=samplename_sorted.bam \
    SORT_ORDER=coordinate \
    CREATE_INDEX=true

To use this BAM file, you also need to create a BAM index, so that software can efficiently access the compressed file. You notice that we automatically include index creation when sorting by specifying `CREATE_INDEX=true`. This is a universal command that can be applied to any PicardTools program. However, if you need to create an index separately, we do this with the **[BuildBamIndex](http://broadinstitute.github.io/picard/command-line-overview.html#BuildBamIndex)** command.

    java -jar ~/path/to/picard.jar BuildBamIndex \
    I=samplename_sorted.bam

Now let's do this for our data. As before, we'll run a test on one file, and then convert to a loop to process all BAMs.

First, we need to load the right module:

    module load jdk/10.0.1-fasrc01
    
We'll also set the $PICARD_HOME variable, so we can just type $PICARD_HOME/picard.jar instead of a long path.

    PICARD_HOME=/n/sw/fasrcsw/apps/Core/picard/2.9.0-fasrc01

Now, let's run our test:

    java -Xmx4g -XX:ParallelGCThreads=1 -jar $PICARD_HOME/picard.jar SortSam \
    I=Falb_COL1.sam \
    O=Falb_COL1.sorted.bam \
    SORT_ORDER=coordinate \
    CREATE_INDEX=true

You can check to make sure you have a bam and a bam index (.bai) file with ls.

Now let's do a loop:

    for INDEX in 1 2 3 4 5 31 32 33 34 35
    do
      java -Xmx4g -XX:ParallelGCThreads=1 -jar $PICARD_HOME/picard.jar SortSam \
      I=Falb_COL$INDEX.sam \
      O=Falb_COL$INDEX.sorted.bam \
      SORT_ORDER=coordinate \
      CREATE_INDEX=true
    done

### Alignment Metrics

It may also be useful to calculate metrics on the aligned sequences. We can easily do this with the **[CollectAlignmentSummaryMetrics](http://broadinstitute.github.io/picard/command-line-overview.html#CollectAlignmentSummaryMetrics)** tool. Note that you can collect the alignment metrics on several different levels. In the below example, I've included metrics both at the sample and read group level. You also need to include the reference fasta file.

    java -jar ~/path/to/picard.jar CollectAlignmentSummaryMetrics \
    I=samplename_sorted.bam \
    R=reference.fasta \
    METRIC_ACCUMULATION_LEVEL=SAMPLE \
    METRIC_ACCUMULATION_LEVEL=READ_GROUP \
    O=samplename.alignment_metrics.txt

In the interest of time, we are not going to run this today on our samples, but you should be able to figure out how to do this yourself from the sample code above. 

### Deduplication

After alignment, sorting and indexing, it is necessary to identify any duplicate sequences from the same DNA fragment in your files that occur due to sample preparation (e.g. during PCR) or incorrect optical cluster identification during sequencing. This possibility is why it is important to identify read groups for different lanes of the same sample. This is also a useful point to merge together any BAM files from the same sample that are currently separated (demonstrated in example below). We identify duplicate sequences with **[MarkDuplicates](https://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates)**, and additional details on how this is performed can be found in the tool documentation. 

Note that it is not recommended to actually remove the duplicate sequences from the file, but simply to mark the flags appropriately in the BAM file, so that those sequences are ignored downstream. If using tools other than those we recommend here, make sure they can identify these flags. These sequences can also be removed later should the need arise.

At this point, we should be familar with the loop command, so we are not going to run a test with just a single sample. Here is the full command:

    for INDEX in 1 2 3 4 5 31 32 33 34 35
    do
      java -Xmx4g -XX:ParallelGCThreads=1 -jar $PICARD_HOME/picard.jar MarkDuplicates \
      TMP_DIR=tmp \
      I=Falb_COL$INDEX.sorted.bam \
      O=Falb_COL$INDEX.dedup.bam \
      METRICS_FILE=Falb_COL$INDEX.dedup.metrics.txt \
      REMOVE_DUPLICATES=false \
      TAGGING_POLICY=All
    done

We also recommend creating a deduplications metrics file, which will report the proportion and type of duplicate sequences in your sample and read groups.

***Following deduplication make sure to sort and index your file again, as shown in the above section.***

    for INDEX in 1 2 3 4 5 31 32 33 34 35
    do
      java -Xmx4g -XX:ParallelGCThreads=1 -jar $PICARD_HOME/picard.jar SortSam \
      I=Falb_COL$INDEX.dedup.bam \
      O=Falb_COL$INDEX.final.bam \
      SORT_ORDER=coordinate \
      CREATE_INDEX=true
    done
jav
### Validating BAM files

Once you are done with the above steps, it is best practice to validate your BAM file, to make sure there were not issues or mistakes associated with previous analyses. This is done with **[ValidateSamFile](http://broadinstitute.github.io/picard/command-line-overview.html#ValidateSamFile)**.

    java -jar ~/path/to/picard.jar ValidateSamFile \
    I=sample.dedup.sorted.bam \
    O=sample.validate.txt \
    MODE=SUMMARY

ls
At this point, we have an analysis-ready BAM. For some workflows, we would do indel realignment and base quality score recalibration at this point. These are described in more detail on the [Informatics Pop Gen Tutorial](https://informatics.fas.harvard.edu/whole-genome-resquencing-for-population-genomics-fastq-to-vcf.html), from which this workshop is derived. 

With larger datasets, you would also want to remove some intermediate files once you have successfully validated your final bams. 

## Variant calling <a name="variantcalling"></a>

There are multiple options for variant calling, including programs like [FreeBayes](https://github.com/ekg/freebayes), [Samtools](http://www.htslib.org), and the [GATK](https://software.broadinstitute.org/gatk/). For this tutorial, we are focusing on the **[HaplotypeCaller](https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_gatk_tools_walkers_haplotypecaller_HaplotypeCaller.php)** program from the GATK pipeline. Calling variants with HaplotypeCaller is essentially a two-step process (similar to indel realignment). First, you call genotypes individually for each sample. Second, you perform joint genotyping across samples to produce a multi-sample VCF call-set. The advantage to this strategy is that the most computationally intensive step, calling genotypes for each sample, only needs to be performed once, even if additional samples will be added later. The joint genotyping, which is less computationally intensive, can be performed as many times as needed as individuals may be added to the dataset.

*Note that even if you are not planning on using SNP calls in downstream analyses (e.g. due to low-coverage sequencing), it is possible to use the genotype likelihood scores from the resulting VCF files (discussed below), and take advantage of the active development on these variant callers.*

### 1. Calling variants for each sample

For each sample, the **[HaplotypeCaller](https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_gatk_tools_walkers_haplotypecaller_HaplotypeCaller.php)** program is used to call variants. The minimum options needed are a reference genome, BAM files for that sample, and output file name. Note that this process can be computationally intensive, so to speed the process up, you may wish to use a  *scatter-gather* approach (Figure 2), and perform the variant calling on non-overlapping segments of the genome, specified with the `-L scaffold_list` [option](https://software.broadinstitute.org/gatk/documentation/article.php?id=4133). 

The output for this program will be a [GVCF](https://software.broadinstitute.org/gatk/documentation/article?id=4017), which has raw, unfiltered SNP and indel calls for all sites, variant or invariant, unlike a typical VCF file (see below for descriptions of the VCF file format.). This is specified by the `--emitRefConfidence GVCF` command, with an example below. See the program page for additional parameter options.

***Note, for low-coverage data, we recommend changing the defaults for two options: `-minPruning 1` and `-minDanglingBranchLength 1`. These commands ensure that any paths in the sample graph (see [detailed documentation on the model](https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_gatk_tools_walkers_haplotypecaller_HaplotypeCaller.php)) are only dropped if there is no coverage. Otherwise the defaults of 2 and 4 respectively, will drop low-coverage regions. See the [documentation](https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_gatk_tools_walkers_haplotypecaller_HaplotypeCaller.php) for details on these and other available options.***

Like Picard, let's make a GATK_HOME variable to simplify typing:

    GATK_HOME=/n/holylfs/LABS/informatics/workshops/read_mapping_variant_calling/gatk-4.0.8.1

For our data, we will first move to the vcfs directory, and we will use the low coverage options:

    cd ../03_vcfs
    for INDEX in 1 2 3 4 5 31 32 33 34 35
    do
     $GATK_HOME/gatk \
     HaplotypeCaller \
     --java-options "-Xmx4g -XX:ParallelGCThreads=1" \
     -R ../00_genome/Falbicolis.chr5.fa \
     -I ../02_bams/Falb_COL$INDEX.final.bam \
     -O Falb_COL$INDEX.raw.g.vcf \
     --emit-ref-confidence GVCF \
     --min-pruning 1 \
     --min-dangling-branch-length 1
    done

This is the step can take a few minutes to run, so let's use this time for some questions if you have them.

### 2. Joint genotyping across samples

Once you have run HaplotypeCaller on your cohort of samples, you can use **[GenotypeGVCFs]()** to perform joint genotyping and produce a multi-sample variant call-set from your gVCF files. 

If you have split up your genome into intervals for Haplotype calling above, you can continue to run variant calling on those intervals by adding the `-L scaffold_list` [option](https://software.broadinstitute.org/gatk/documentation/article.php?id=4133). Also, note that for non-human organisms, you may wish to vary the heterozygosity prior from the default value of *0.001*. You can do this with the heterozygosity option, for example with a  value of *0.005*: `--heterozygosity 0.005`. If you wish to include all sites, both variant and invariant, you need to use the `--includeNonVariantSites true` option. See the program page for additional parameter options.

Here is what we'll run. Note that this is not a loop, as we are going to make just one final VCF for all samples.


######## GVCF TO COMBINED VCF CODE

$GATK_HOME/gatk \
GenomicsDBImport \
--java-options "-Xmx4g -XX:ParallelGCThreads=1" \
-R ../00_genome/Falbicolis.chr5.fa \
-V Falb_COL1.raw.g.vcf \
-V Falb_COL2.raw.g.vcf \
-V Falb_COL3.raw.g.vcf \
-V Falb_COL4.raw.g.vcf \
-V Falb_COL5.raw.g.vcf \
-V Falb_COL31.raw.g.vcf \
-V Falb_COL32.raw.g.vcf \
-V Falb_COL33.raw.g.vcf \
-V Falb_COL34.raw.g.vcf \
-V Falb_COL35.raw.g.vcf \
--genomicsdb-workspace-path my_database \
--intervals 5:21000000-23000000

$GATK_HOME/gatk \
GenotypeGVCFs \
--java-options "-Xmx4g -XX:ParallelGCThreads=1" \
-R ../00_genome/Falbicolis.chr5.fa \
-V gendb://my_database \
-O Falb.final.vcf \
--heterozygosity 0.005

*Note, if you specify output file names with `.gz` extensions, GATK will automatically compress your output files and create and index with tabix.*

### VCF File Format

A standard way to store variant information in single or multiple samples is a **[Variant Call Format, or VCF file](https://software.broadinstitute.org/gatk/documentation/article?id=1268)**. The general format is a header, with information about the dataset, references, and annotations, and these lines start with `##`. Following the header, each variant is represented on a separate line with tabs separating each field. It starts with information on the chromosome (CHROM), position (POS), variantID (ID), reference allele (REF), alternate allele (ALT), quality score (QUAL), filter designation (FILTER, e.g. PASS), annotations (INFO),  individual representation format (FORMAT), and finally genotype information for each sample.

The sample-level information can vary depending on the program used, but with the GATK pipeline it is shown as: `GT:AD:DP:GQ:PL`.
Where: 

* `GT` : **Genotype**<br>
    Genotype for the sample at each site. For a diploid, 0/0 indicates homozygous reference alleles, 0/1 indicates a heterozygous sample, with a reference and alternate allele, 1/1 indicates homozygous alternate allele, and ./. indicates no sequence at that site for that individual. Note that samples with different ploidies will have an appropriate number of alleles.

* `AD`: **Allele Depth**<br>
    The number of reads that support each allele. If diploid, ref,alt.

* `DP`: **Depth of Coverage**<br>
    The filtered depth of coverage at the sample level.
    
* `GQ`: **Genotype Quality**<br>
    The quality score of the genotype for that individual. This is usually measured on the [Phred scale](https://software.broadinstitute.org/gatk/documentation/article.php?id=4260), as with the Fastq file format, described above. Higher values are better.

* `PL`: **"Normalized" Phred-scaled likelihoods of possible genotypes**<br>
    For each possible genotype (three in the case of a diploid), the normalized likelihood score ([phred](https://software.broadinstitute.org/gatk/documentation/article.php?id=4260)-scaled), with the most likely genotype set to 0. The other values are scaled to this most likely genotype.

## Data filtering <a name="filtering"></a>

Once you have produced a VCF file with all of your samples, it is necessary to [evaluate](https://gatkforums.broadinstitute.org/gatk/discussion/6308/evaluating-the-quality-of-a-variant-callset) and filter the dataset for quality. If you are working with an organism that has a **variant truth set**, or set variants that are thought to be correct, prior to any downstream analyses you should perform **[Variant Quality Score Recalibration (VQSR)](https://gatkforums.broadinstitute.org/gatk/discussion/39/variant-quality-score-recalibration-vqsr)**. As VQSR is a tricky process to get right (and with few organisms outside of humans with appropriate truth and training datasets), the GATK has a [detailed tutorial](https://software.broadinstitute.org/gatk/documentation/article.php?id=2805). Here, we will focus on using **[hard filters](https://software.broadinstitute.org/gatk/documentation/article.php?id=3225)**, with a tutorial based on the [GATK recommendations](https://software.broadinstitute.org/gatk/documentation/article.php?id=2806).

We will not have time to get into filtering today, but you can read more about filtering options in the [Informatics Pop Gen Tutorial](https://informatics.fas.harvard.edu/whole-genome-resquencing-for-population-genomics-fastq-to-vcf.html).
