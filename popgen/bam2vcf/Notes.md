Preliminary Steps 
1. Connect to VPN
2. Connect to virtual desktop
	- type: 'vdi.rc.fas.harvard.edu' into address bar in browser
	- or google FAS Harvard vdi and follow links
3. Click 'Interactive Apps' tab, and select Jupyter Lab
4. Wait for connection, then press blue 'Connect' box
5. Create a Terminal tab (and type "bash" so your prompt shows the current directory you're in) and a text file tab
	- we will use the text file to take notes and write commands
 	- we will copy these commands into the terminal window



Today we will be talking about calling variants, which include SNPs and indels. Pipelines like the one we'll walk through today need to be used any time you're interested in detecting mutational differences between samples, compared to using sequence reads to calculate relative abundance (e.g. RNAseq differential expression, CHIPseq, ATACseq)

For the sake of efficiency we will work in the same environment that previous sessions of this workshop used in which various programs are already installed. This will give us more time to explore the functionality of these programs, but unfortunately skips the step of installing or setting up these programs in the first place. I've ran this workflow on the cluster (as opposed to within this virtual environment) and will try to explain how I set up various programs.

## Programs used in this session
### Modules: samtools, GATK4 (+java8), vcftools
### on odyssey: Picard 
### NOT on odyssey: snpEff

The first 3 programs are available as modules on odyssey. For instance, if we Google "Harvard FAS module", we can search for the commands to load many of these tools. After they are loaded, you may simple run these programs by typing the name.


However, since the GATK is updated so frequently, many may choose to download the newest version directly from the website and place it in your home directory. For instance, we can Google GATK and go to the download tab on top of the homepage. To run GATK, just put the entire path/location (where you have put the program) before the program name, e.g. "/n/home/bjarnold/gatk HaplotypeCaller ...".

When using GATK on the cluster, it's important to also use a specific version of Java that it's compatible with, version 8, not version 9. Simpy type 'module load jdk/1.8.0_45-fasrc01' to load a version of java that is compatible with GATK.


To annotate VCFs we will use SnpEff, which currently is not available on odyssey, but again you can download this program from the website and it has instructions on how to set it up and verify that it's working. Again, if you install this yourself, you can run the program by specifying the full file path, like 
```
java -Xmx4g -jar /n/home11/bjarnold/snpEff/snpEff.jar -v loxAfr3.86 sample.vcf > sample_snpEff.vcf.
```
where we type "java" to run this program, just as we type "python" before python scripts we want to run.


Let's get started! to run GATK, we need two things: the reference genome to which we mapped reads, and the mapped reads in BAM file format. You should have both of these from the last workshop session, but I provide them here in cased you missed it.

```
DATA_DIR=/n/home11/bjarnold/workshops/bam2vcf_2019
ls $DATA_DIR
```

where the bam files are in the "bams" directory. there's also a file called "PicardPath.txt" that contains the file path to picard, which is on odyssey but not loaded as a module. Lets open this file with "cat" and copy-paste this path into our notebook

```
cat $DATA_DIR/picardPath.txt
```


### Process reference genome
Even though we have our reference genome, the GATK requires processing your reference genome to make 2 additional files that allow it to efficiently access the bases. Specifically, we need to create an index and a dictionary.


```
mkdir gatkWorkshop2019
cd gatkWorkshop2019
mkdir 00_genome
cp $DATA_DIR/Falbicolis.chr5.fa 00_genome
cd 00_genome
PICARD_HOME=/n/sw/fasrcsw/apps/Core/picard/2.9.0-fasrc01
java -jar $PICARD_HOME/picard.jar CreateSequenceDictionary \
R=Falbicolis.chr5.fa \
O=Falbicolis.chr5.dict
samtools faidx Falbicolis.chr5.fa
cd ..
```

Make sure you have 3 files! a .fa (reference), .fai (the index), and .dict (the dictionary) file

## VARIANT CALING WITH GATK's HAPLOTYPE CALLER
### GATK HaplotypeCaller workflow: BAM -> GVCF -> DB -> VCF

HaplotypeCaller (HC) calls SNPs and indels via local de-novo reassembly of haplotypes. When the program encounters a region of polymorphism, it throws away the existing mapping variation and reassembles the reads in the region to perform assembly and variant calling simultaneously. This allows for more accurate variant detection in regions that are traditionally difficult to call, e.g. when multiple variants (SNPs+indels) are very close to one another.

HC is not the only variant caller, and all methods have their strengths and weaknesses. Compared to others, HC is more difficult to use. However, while a program like Samtools may perform quite well for SNP calling, the ability for HC to simultaneously call indels is a fantastic feature. When I was testing out various programs that specialize in indel calling (e.g. svABA), HaplotypeCaller gave extremely similar results, although this depends on the size of the indels you want to detect.

The GATK HC pipeline contains multiple steps designed to streamline large sequencing experiments in which new data continually comes in. An intermediate GVCF file is first created per sample, and these individual GCVFs are then combined across samples. If more samples are sequenced, the previous samples need not go thorugh this initial step again.

After the first step of creating per-sample GVCFs, two additional steps are used to combine this information by creating a database that is then used to jointly genotype the entire sample.

### HAPLOTYPE CALLER: BAM -> GVCF
This produces our intermediate GVCF files which contains information about ALL sites, whether or not there is a variant there, which is important for when joint genotyping is performed across all samples. I typically never look at these files. This step may produce a lot of WARNINGS, which are different from ERRORS and are usually completely normal. There are ways to silence these warnings, but we won't do this.

```
mkdir 01_gvcfs
cd 01_gvcfs
```

Lets make some variables to simplyfy typing and store directory names

```
BAM_DIR=$DATA_DIR/bams
REF_FILE=../00_genome/Falbicolis.chr5.fa
```

We can verify our BAM files are within this directory, and instead of typing out the gatk command for each one, we can create a loop. Here we will create a loop that goes across all samples, but it may be useful to only run one sample at first just to make sure it works!

```
for INDEX in 1 2 3 34 35
do
	gatk HaplotypeCaller \
	--java-options "-Xmx4g -XX:ParallelGCThreads=1" \
	-R $REF_FILE \
	-I $BAM_DIR/Falb_COL$INDEX.final.bam \
	-O Falb_COL$INDEX.raw.g.vcf \
	-L 5 \
	--emit-ref-confidence GVCF \
	--min-pruning 1 \
	--min-dangling-branch-length 1 
done
```
the interval we specified corresponds to the only scaffold in the fasta file, named "5"
```
head $REF_FILE
```

Note: for low-coverage data, we recommend changing the defaults for two options: -minPruning 1 and -minDanglingBranchLength 1. These commands ensure that any paths in the sample graph (see detailed documentation on the model) are only dropped if there is no coverage. Otherwise the defaults of 2 and 4 respectively, will drop low-coverage regions. See the documentation for details on these and other available options.

The -L argument lets GATK know that you only want to condiser a specific interval, here chromosome 5, which corresponds to the name of the only sequence in our reference. You could make another loop that iterates over chromosome numbers, but this depends on how chromosomes or scaffolds are named. Alternatively, you can create a file that contains a list of intervals, and just give that to GATK instead. If you're dealing with bacterial genomes that are relatively small, you probably don't need this, and not specifying -L will just do the whole genome. However, for larger eukaryotic genomes you'll certainly need to do this if you have lots of data and don't want this program to take weeks to run.

### GenomicsDBImport: GVCF -> DB
Before we can jointly genotype the GVCFs we just produced, we must first merge these GVCFs into a genomics database that stores information in a specialized way. For more information about this, see the GATK website, but the purpose of this preliminary gathering of GVCFs is to enhance computational speed. 

GenomicsDBImport uses temporary disk storage during import. The amount of temporary disk storage required can exceed the space available, so we can specify a specific 'tmp' directory that we know has sufficient space, e.g. on scratchlfs.

GenomicsDBImport also requires specifying an interval, which could be typed in directly or you could have a file in which all of the intervals you want are on separate lines.

```
cd ..
mkdir -p 02_genDB/tmp
cd 02_genDB
GVCF_DIR=../01_gvcfs

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

### GenotypeGVCFs: DB -> VCF
Finally, lets generate our VCF file.

```
cd ..
mkdir -p 03_vcfs/tmp
cd 03_vcfs
DB_DIR=../02_genDB

gatk GenotypeGVCFs \
--java-options "-Xmx4g" \
-R $REF_FILE \
-V gendb://$DB_DIR/my_database \
-O final.vcf \
--tmp-dir=tmp
```

### THE VCF FORMAT

A standard way to store variant information is a Variant Call Format, or VCF file. The general format is a header, 
#with information about the dataset, references, and annotations, and these lines start with ##. The last line of the header, which starts with '#CHROM', shows how variant information is displayed.

```
grep '#CHROM' final.vcf
```

Following the '##' header, the line with '#CHROM' shows how information for each variant is represented on a separate line with tabs separating each field. It starts with information on the chromosome (CHROM), position (POS), variantID (ID), reference allele (REF), alternate allele (ALT), quality score (QUAL), filter designation (FILTER, e.g. PASS), annotations (INFO), individual representation format (FORMAT), and finally genotype information for each sample.

To get an idea of what these columns represent, let's look a tthe first two variants in our VCF, the two that come right after the header. We could just open the VCF file in a text editor, but oftentimes these files are huge, making this not practical. But, but we can easily glance at particular lines using UNIX commands.

```
grep -v '#' final.vcf | head -n 2
```

The INFO column contains a lot of information about the data that went into calling this variant. We'll discuss this shortly in the next section on filtering.Unlike the INFO field that has information for the entire site, the FORMAT column shows how information for each individual is represented in the genotype column, it is shown as: GT:AD:DP:GQ:PL. Where:

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
We currently have a VCF file, but it would also be useful to know, for each variant, where it occurs in the genome. Is it within a gene or is it intergenic? If it's within a gene, does it change the amino acid sequence (nonsynonymous/missense)?

We can use the program snpEff to add this information to our VCF file. Let's use another VCF file that I provide, as opposed to the one we just made for the Collared Flycatcher birds. This new VCF is just the first 50k lines of a much larger file on savanah elephants. It currently is not annotated, but we can add this information using the following snpEff command

```
cd ..
mkdir 04_snpeff
cd 04_snpeff
cp $DATA_DIR/subsampled.vcf ./
snpEff -dataDir $PWD -v loxAfr3.86 subsampled.vcf > subsampled_snpEff.vcf
```

where the -v command indicated "verbose" mode, which I like because this information can help debug anything. After you list all your options with dashes, you need to type the name of the database you want to use, followed by the name of the VCF file we want to annotate. Then, we store this output into a file that will be our annotated VCF. Here I use the Loxodonta Africana database, which is the genome I aligned my reads to, and you can easily search which databases are available within snpEff by either accessing the .config file or by typing "snpEff databases", which outputs a lot of text so you may want to direct that to a file. You can use this to find the ID of your database, e.g. loxAfr3.86.

### BASIC FILTERING/SANITY CHECKS WITH VCFTOOLS
For many basic tasks that involve VCF files, chances are high that a program has already been developed to do just that. Vcftools is a widely-used program that has numerous methods for working with VCF files, including ways to further filter the VCF file, e.g. according to allele frequencies, sequencing depths, etc., and also provides tools to convert your VCF file to other popular file formats such as PLINK. You should certianly read the manual for this program, but we will use several of its features to demonstrate its utility.

For instance, one common task is to see if there are particular samples in the data that are just bad, and have a lot of missing data perhaps from lower sequencing depths or from contamination, which causes fewer reads to map to the reference genome provided. It would be good to identify these individuals and exclude them from downstream analyses for the following reason. Let's say you want to do an analysis that uses genotypes, and you want only sites in which each individual in the sample has at least 4X depth to ensure genotypes are called with sufficient information, and any site that has at least one individual with less depth than this threshold gets excluded. Thus, if one individual has dramatically less sequencing depth than all others, it may cause many sites to get filtered out, such that it may be best to instead exclude that individual from the analysis altogether.

We can get a sense of the amount of missing data per individual using the following options in Vcftools:
```
vcftools --vcf subsampled_snpEff.vcf --remove-filtered-all --minDP 4 --out subsampled --missing-indv
```

Looking at subsampled.imiss, T1A has a ton of missing data... turned out this individual was highly contaminated because tissue was sampled from a dead elephant.

Some analyses, especially in poplation genetics, assume you have a random sample from a population, but what if we unknowingly sampled individuals witin a family. We can check this using the 'relatedness' option with Vcftools.

```
vcftools --vcf subsampled_snpEff.vcf --remove-filtered-all --minDP 4 --remove-indv T1A --out subsampled --relatedness
```

A relatedness score of 1 is expected when comparing an individual with itself, and a score of 0 is expected between individuals randomly sampled from a population. These values all look a little high and strange, but this is because we're working with so few sites. Nonetheless, one particular pair of individuals does look to be quite related, and this is confirmed when analyzing the entire VCF with many millions of variants.

Lastly, another popular program to analyze sequence data is PLINK. I will not cover PLINK here, but I will note that it has a specialized input file format. We can convert our VCF file into PLINK format using vcftools
```
vcftools --vcf subsampled_snpEff.vcf --remove-filtered-all --minDP 4 --remove-indv T1A --out subsampled --plink
```


### BRIEF INTRODUCTION TO CUSTOM VCF ANALYSES WITH CYVCF2
For much of the basic analyses you'd like to do, existing software will get you most of the way there. But what if there's something else you'd like to do with your VCF file.

Because we know the structure of a VCF file, we can easily write a python script to get the information we want on each line by just splitting based on tabs, then we can further split on colons or pipes based on what information we want in individual fields. However, a python module does just this, which will make our code look much more simple and save us some time coding. Several python modules exist to parse python modules, and I would recommend *ONLY* using cyvcf2; the others that I have used are frustratingly slow and can take days or weeks to get through large VCF files. Cyvcf2 is extremely fast!!

We've pre-installed this module for today so that we don't have to wait 15 minutes for it to install, but you can easily install it using conda, and I suggest you make a specific conda environment for it.

Once that's done, here is a very simple python script you may use as a template for your own custom analyses.

# *check out simple_cyvcf2_ex.py*

