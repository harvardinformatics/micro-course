Preliminary Steps 
1. Connect to VPN
2. Connect to virtual desktop
	- type: 'vdi.rc.fas.harvard.edu' into address bar in browser
	- or google FAS Harvard vdi and follow links
3. Click 'Interactive Apps' tab, and select Jupyter Lab
4. Wait for connection, then press blue 'Connect' box
5. Create a Terminal tab (and type "bash" so your prompt shows the current directory you are in) and a text file tab
	- we will use the text file to take notes and write commands
 	- we will copy these commands into the terminal window



This workshop is about mapping short sequence reads to a genome, and using that data to call variants, which include SNPs and small indels. Pipelines like the one we will walk through today need to be used any time you arre interested in detecting mutational differences between samples, compared to using sequence reads to calculate relative abundance (e.g. RNAseq differential expression, CHIPseq, ATACseq).

We will work in the same environment that previous sessions of this workshop used in which various programs are already installed. This will give us more time to explore the functionality of these programs, as everything is already set up and working, but unfortunately skips the step of installing or setting up these programs in the first place. I have ran this workflow on the cluster (as opposed to within this virtual environment) and will try to explain how I set up various programs.

## Programs used in this session
### Modules: samtools, bwa, GATK4 (+java8 !!), vcftools
### on odyssey: Picard 
### NOT on odyssey: snpEff

The first 3 programs are available as modules on odyssey. For instance, if we Google "Harvard FAS module", we can search for the commands to load many of these tools. After they are loaded, you may simply run these programs by typing the name.

If you remember that these modules are on Cannon but forgot exactly what to type, you can always type 'module-query' followed by the name of the program, e.g. samtools, to find exactly how to load the module:

```
module-query samtools
```

For programs that are not available as modules, we can download them and place them in our home directory. For instance, we can Google GATK and go to the download tab on top of the homepage to get the newest version, which may not be on Cannon. To then run GATK, just put the entire path/location (where you have put the program) before the program name, e.g. "/n/home/bjarnold/gatk HaplotypeCaller ...".

Sometimes, to run programs you have to type a little more than just the name of the program. For instance, to annotate VCFs we will use SnpEff, which currently is not available on Cannon, but again you can download this program from the website and it has instructions on how to set it up and verify that it is working. Again, if you install this yourself, you can run the program by specifying the full file path, like 

```
java -Xmx4g -jar /n/home11/bjarnold/snpEff/snpEff.jar -v loxAfr3.86 sample.vcf > sample_snpEff.vcf.
```

where we type "java" to run this program, just as we type "python" before python scripts we want to run. Running other programs represented as .jar files will follow a similar syntax, like picard.


Let's get started! All of the data we'll use is in the following directory:

```
DATA_DIR=/n/home11/bjarnold/workshops/bam2vcf_2019
ls $DATA_DIR
```

Here we have two directories that we will use to start the workshop that contain the reference genome, and the raw Illumina reads that we will map to this genome. There is also a file called "PicardPath.txt" that contains the file path to picard, which is on odyssey but not loaded as a module. Lets open this file with "cat".

```
cat $DATA_DIR/picardPath.txt
```

Today we will not need to do this because we have specially installed picard for this workshop, but normally to get picard working on the cluster you have to run the java .jar file with

```
java -Xmx4g -jar $PICARD_HOME/picard.jar SortSam
```

Where SortSam is the specific Picard command you may want to run.




Lets make a directory for all the work we will do for this workshop and copy the relevant files into it.

```
mkdir readMapVarCall
cd readMapVarCall
cp -r $DATA_DIR/00_genome ./
cp -r $DATA_DIR/01_fastqs ./
```

## ALIGNING READS TO REFERENCE GENOME WITH BWA
### Fastq -> SAM -> BAM -> Mark duplicates -> INDEX

### Process reference genome
Both bwa (and GATK) requires processing your reference genome to make additional files that allow it to efficiently access the bases.

```
cd 00_genome
bwa index -p Falb Falbicolis.chr5.fa
cd ..
```

You should have several new files that all start with 'Falb', since this is the value we gave after the '-p' flag.

Lets use BWA to align our reads to the reference genome and get a SAM file (Sequence Alignment/Map), one for each sample. For information about BWA, simply type its name, followed by "--help". This procedure works for most programs.

```
mkdir 02_bams
cd 02_bams
for INDEX in 1 2 3 34 35;
do
   bwa mem -M -t 1 -R "@RG\tID:COL_$INDEX\tSM:COL_$INDEX" ../00_genome/Falb \
   ../01_fastqs/Falb_COL$INDEX.1.fastq \
   ../01_fastqs/Falb_COL$INDEX.2.fastq \
   > Falb_COL$INDEX.sam 2> Falb_COL$INDEX.log
done
```

Note: read groups refer to sets of reads that were generated from a single run of a sequencing instrument, and this information is stored in our SAM file on lines that start with "@RG". Using read groups allows us to not just distinguish between samples, but also particular samples that were sequenced across several experiments. Programs like the GATK require this information so that it can attempt to compensate for variability between sequencing runs.

Next, lets convert our SAM files from BWA to BAM files, which are compressed versions that a lot of downstream programs use as input files.

```
for INDEX in 1 2 3 34 35;
do
  picard SortSam \
  I=Falb_COL$INDEX.sam \
  O=Falb_COL$INDEX.sorted.bam \
  SORT_ORDER=coordinate \
  CREATE_INDEX=true
done
```

In preparation for the GATK, we need to mark paired-end reads that represent PCR duplicates. For instance, if you randomly sheared your DNA before sequencing, the chance the both the forward and reverse read have the exact same coordinates as another pair of reads is highly unlikley. A more likely explanation is that the same molecule in your DNA library was duplicated during the PCR step. This can potentially be a problem if there were particular molecules that, by chance, replicated many times. In any case, these duplicates dont represent random samples of DNA fragments are are usually flagged prior to running GATK. 

However, if you are using sequencing protocols in which you expect many read pairs to have the same start and stop positions (e.g. RADseq), do not do this step!!

```
for INDEX in 1 2 3 34 35
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

Last, we need to create an index of our BAM file in order for downstream programs to quickly access its contents.

```
for INDEX in 1 2 3 34 35
do
  picard BuildBamIndex \
  I=Falb_COL$INDEX.dedup.bam
done
cd ..
```

## VARIANT CALING WITH GATK HAPLOTYPE CALLER
### GATK HaplotypeCaller workflow: BAM -> GVCF -> DB -> VCF

For GATK, we also need to process the reference genome. Specifically, we need to create an index and a dictionary.

```
cd 00_genome
picard CreateSequenceDictionary \
R=Falbicolis.chr5.fa \
O=Falbicolis.chr5.dict
samtools faidx Falbicolis.chr5.fa
cd ..
```

Make sure you have a .fai (the index) and .dict (the dictionary) file.

Some background:

HaplotypeCaller (HC) calls SNPs and small indels via local de-novo reassembly of haplotypes. When the program encounters a region of polymorphism, it throws away the existing mapping variation and reassembles the reads in the region to perform assembly and variant calling simultaneously. This allows for more accurate variant detection in regions that are traditionally difficult to call, e.g. when multiple variants (SNPs+indels) are very close to one another.

HC is not the only variant caller, and all methods have their strengths and weaknesses. Compared to others, HC is more difficult to use. However, while a program like Samtools may perform quite well for SNP calling, the ability for HC to simultaneously call indels is a fantastic feature. However, other programs are starting to have this capability, but GATK is nonetheless an extremely popular variant caller.

The GATK HC pipeline contains multiple steps designed to streamline large sequencing experiments in which new data continually comes in. An intermediate GVCF file is first created per sample, and these individual GCVFs are then combined across samples. If more samples are sequenced, the previous samples need not go thorugh this initial step again.

After the first step of creating per-sample GVCFs, two additional steps are used to combine this information by creating a database that is then used to jointly genotype the entire sample.

### Part 1: HAPLOTYPE CALLER: BAM -> GVCF
This produces our intermediate GVCF files which contains information about ALL sites, whether or not there is a variant there, which is important for when joint genotyping is performed across all samples. I typically never look at these files. This step may produce a lot of WARNINGS, which are different from ERRORS and are usually completely normal. There are ways to silence these warnings, but we will not do this.

```
mkdir 03_gvcfs
cd 03_gvcfs
```

Lets make some environmental variables to simplify typing and store directory names:

```
BAM_DIR=../02_bams
REF_FILE=../00_genome/Falbicolis.chr5.fa
```

If these environmental variables are not specified, they contain no values and GATK will complain that youve specified an invalid argument.

We can verify our BAM files are within this directory, and that we have correctly specified the path to the reference file. 

```
ls $BAM_DIR
head $REF_FILE
```

Instead of typing out the gatk command for each BAM file, we can create a loop that does this for us. Here we will create a loop that goes across all samples, but it may be useful to only run one sample at first just to make sure it works!

```
for INDEX in 1 2 3 34 35
do
	gatk HaplotypeCaller \
	--java-options "-Xmx4g -XX:ParallelGCThreads=1" \
	-R $REF_FILE \
	-I $BAM_DIR/Falb_COL$INDEX.dedup.bam \
	-O Falb_COL$INDEX.raw.g.vcf \
	-L 5 \
	--emit-ref-confidence GVCF \
	--min-pruning 1 \
	--min-dangling-branch-length 1 
done
```

The interval we specified , via -L, corresponds to the only scaffold in the fasta file, named "5".

Note: for low-coverage data, we recommend changing the defaults for two options: -minPruning 1 and -minDanglingBranchLength 1. These commands ensure that GATK still considers sites with low sequencing depth, otherwise it drops these regions. The defaults are 2 and 4 respectively. See the documentation for details on these and other available options.

The -L argument lets GATK know that you only want to condiser a specific interval, here chromosome 5, which corresponds to the name of the only sequence in our reference. You could make another loop that iterates over chromosome numbers, but this depends on how chromosomes or scaffolds are named. Alternatively, you can create a file that contains a list of intervals, and just give that to GATK instead. If you are dealing with bacterial genomes that are relatively small, you probably do not need this, and not specifying -L will just do the whole genome. However, for larger eukaryotic genomes you will certainly need to do this if you have lots of data and do not want this program to take weeks to run.

### Part 2: GenomicsDBImport: GVCF -> DB
Before we can jointly genotype the GVCFs we just produced, we must first merge these GVCFs into a genomics database that stores information in a specialized way. For more information about this, see the GATK website, but the purpose of this preliminary gathering of GVCFs is to enhance computational speed for combining genotype information across samples. 

GenomicsDBImport uses temporary disk storage during import. The amount of temporary disk storage required can exceed the space available, so we can specify a specific 'tmp' directory that we know has sufficient space, e.g. on scratchlfs. I say this because specifying a 'tmp' directory solved a memory issue I once had.

GenomicsDBImport also requires specifying an interval, which could be typed in directly or you could have a file in which all of the intervals you want to include are on separate lines.

```
cd ..
mkdir -p 04_genDB/tmp
cd 04_genDB
GVCF_DIR=../03_gvcfs

gatk GenomicsDBImport \
--java-options "-Xmx4g" \
--genomicsdb-workspace-path my_database \
-L 5 \
--tmp-dir=tmp \
-V $GVCF_DIR/Falb_COL1.raw.g.vcf \
-V $GVCF_DIR/Falb_COL2.raw.g.vcf \
-V $GVCF_DIR/Falb_COL3.raw.g.vcf \
-V $GVCF_DIR/Falb_COL34.raw.g.vcf \
-V $GVCF_DIR/Falb_COL35.raw.g.vcf 
```

### Part 3:  GenotypeGVCFs: DB -> VCF
Finally, lets generate our VCF file.

```
cd ..
mkdir -p 05_vcfs/tmp
cd 05_vcfs
DB_DIR=../04_genDB

gatk GenotypeGVCFs \
--java-options "-Xmx4g" \
-R $REF_FILE \
-V gendb://$DB_DIR/my_database \
-O final.vcf \
--tmp-dir=tmp
```

### THE VCF FORMAT

A standard way to store variant information is a Variant Call Format, or VCF file. The general format is a header, with information about the dataset, references, and annotations, and these lines start with ##. The last line of the header, which starts with '#CHROM', shows how variant information is displayed.

```
grep '#CHROM' final.vcf
```

Following the '##' header, the line with '#CHROM' shows how information for each variant is represented on a separate line with tabs separating each field. It starts with information on the chromosome (CHROM), position (POS), variantID (ID), reference allele (REF), alternate allele (ALT), quality score (QUAL), filter designation (FILTER, e.g. PASS), annotations (INFO), individual representation format (FORMAT), and finally genotype information for each sample.

To get an idea of what these columns represent, we will look a the first two variants in our VCF, the two that come right after the header. We could just open the VCF file in a text editor, but oftentimes these files are huge, making this not practical. But, but we can easily glance at particular lines using UNIX commands.

```
grep -v '#' final.vcf | head -n 2
```

The INFO column contains a lot of information about the data that went into calling this variant. We will discuss this shortly in the next section on filtering. Unlike the INFO field that has information for the entire site, the FORMAT column shows how information for each individual is represented in the genotype column, it is shown as: GT:AD:DP:GQ:PL. Where:

```
GT : Genotype
Genotype for the sample at each site. For a diploid, 0/0 indicates homozygous reference alleles, 0/1 indicates a heterozygous sample, with a reference and alternate allele, 1/1 indicates homozygous alternate allele, and ./. indicates no sequence at that site for that individual. 

AD: Allele Depth
The number of reads that support each allele. If diploid, ref,alt.

DP: Depth of Coverage
The filtered depth of coverage at the sample level.

GQ: Genotype Quality
The quality score of the genotype for that individual. This is usually measured on the Phred scale, as with the Fastq file format, described above. Higher values are better.

PL: stands for "Phred-scaled Likelihoods", which are “Normalized” probabilities of the possible genotypes. For each possible genotype (three in the case of a diploid), the normalized likelihood score (phred-scaled), with the most likely genotype set to 0. The other values are scaled to this most likely genotype.
```

### VCF HARD FILTERING: updating FILTER field using the INFO field
VCF files have a column entitled "FILTER", which at the moment is populated with an uninformative '.' for each site. Within the VCF INFO field is a ton of information about this site, such as the mapping quality of reads that overlap this site (MQ), along with other statistics that may indicate a variant is supported by suspicious evidence, for instance if a variant is only supported by the end of sequencing reads which are known to be enriched for sequencing errors (ReadPosRankSum), or if a variant is supported more by the forward DNA strand but not the reverse strand (FS).

We can use this FILTER column and put a 'PASS' label for any variant that has sufficiently good quality metrics, and the GATK website has a number of recommendations for which metrics to use for filtering. Please see the informatics webpage on variant calling for a link to this information, or go directly to the website by googling something like GATK 'hard filtering'.

While there are more metrics one may use for filtering, here we use only a few just to illustrate how one might hard filter a VCF.

```
gatk VariantFiltration \
-R $REF_FILE \
-V final.vcf \
-O final_filtered.vcf \
--filter-expression "QD < 2.0 || MQ < 40.0 || ReadPosRankSum < -8.0 || FS > 60.0" \
--filter-name "bad"
```

Exercise: using unix commands, how would we get an idea of how many sites are PASS and how many sites were filtered out? Use grep, pipe, wc -l




### ANNOTATION WITH SNPEFF
We currently have a VCF file, but it would also be useful to know, for each variant, where it occurs in the genome. Is it within a gene or is it intergenic? If it is within a gene, does it change the amino acid sequence (nonsynonymous/missense)?

We can use the program snpEff to add this information to our VCF file. We will use another VCF file that I provide, as opposed to the one we just made for the Collared Flycatcher birds. This new VCF is just the first 50k lines of a much larger file on savanah elephants. It currently is not annotated, but we can add this information using the following snpEff command

```
cd ..
mkdir 04_snpeff
cd 04_snpeff
cp $DATA_DIR/subsampled.vcf ./
snpEff -dataDir $PWD -v loxAfr3.86 subsampled.vcf > subsampled_snpEff.vcf
```

where the -v command indicated "verbose" mode, which I like because this information can help debug anything. After you list all your options with dashes, you need to type the name of the database you want to use, followed by the name of the VCF file we want to annotate. Then, we store this output into a file that will be our annotated VCF. Here I use the Loxodonta Africana database, which is the genome I aligned my reads to, and you can easily search which databases are available within snpEff by either accessing the .config file or by typing "snpEff databases", which outputs a lot of text so you may want to direct that to a file. You can use this to find the ID of your database, e.g. loxAfr3.86.

Last, lets convert our VCF file to a table that is much more convenient to work with.

```
ELE_REF=$DATA_DIR/Loxodonta_africana.loxAfr3.dna.toplevel.fa

gatk VariantsToTable \
-R $ELE_REF \
-V subsampled_snpEff.vcf \
-F CHROM -F POS -F QUAL -F AC -F HOM-REF -F HET -F HOM-VAR -F ANN \
-O results.table
```

While many programs you use may accept VCF files as input format, sometimes you may need to do custom analyses with your data. In these cases, it is MUCH easier to make a simple script/program to parse this allele table we created, as opposed to writing a script that parses VCF files which have a much more complicated format.


### BASIC FILTERING/SANITY CHECKS WITH VCFTOOLS
For many basic tasks that involve VCF files, chances are high that a program has already been developed to do just that. Vcftools is a widely-used program that has numerous methods for working with VCF files, including ways to further filter the VCF file, e.g. according to allele frequencies, sequencing depths, etc., and also provides tools to convert your VCF file to other popular file formats such as PLINK. You should certianly read the manual for this program, but we will use several of its features to demonstrate its utility.

For instance, one common task is to see if there are particular samples in the data that are just bad, and have a lot of missing data perhaps from lower sequencing depths or from contamination, which causes fewer reads to map to the reference genome provided. It would be good to identify these individuals and exclude them from downstream analyses for the following reason. For instance say you want to do an analysis that uses genotypes, and you want only sites in which each individual in the sample has at least 4X depth to ensure genotypes are called with sufficient information, and any site that has at least one individual with less depth than this threshold gets excluded. Thus, if one individual has dramatically less sequencing depth than all others, it may cause many sites to get filtered out, such that it may be best to instead exclude that individual from the analysis altogether.

We can get a sense of the amount of missing data per individual using the following options in Vcftools:

```
vcftools --vcf subsampled_snpEff.vcf --remove-filtered-all --minDP 4 --out subsampled --missing-indv
```

Looking at subsampled.imiss, T1A has a ton of missing data... it turned out this individual was highly contaminated because tissue was sampled from a dead elephant.

Some analyses, especially in poplation genetics, assume you have a random sample from a population, but what if we unknowingly sampled individuals within a family. We can check this using the 'relatedness' option with Vcftools.

```
vcftools --vcf subsampled_snpEff.vcf --remove-filtered-all --minDP 4 --remove-indv T1A --out subsampled --relatedness
```

A relatedness score of 1 is expected when comparing an individual with itself, and a score of 0 is expected between individuals randomly sampled from a population. These values all look a little high and strange, but this is because we're working with so few sites. Nonetheless, one particular pair of individuals does look to be quite related, and this is confirmed when analyzing the entire VCF with many millions of variants.

Lastly, another popular program to analyze sequence data is PLINK. I will not cover PLINK here, but I will note that it has a specialized input file format. We can convert our VCF file into PLINK format using vcftools

```
vcftools --vcf subsampled_snpEff.vcf --remove-filtered-all --minDP 4 --remove-indv T1A --out subsampled --plink
```
