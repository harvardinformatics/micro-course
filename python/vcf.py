#!/usr/bin/python -tt

"""
This script goes through a filtered VCF file that has been annotated with snpEff and, for a specified set of groups, calculates differences in allele frequency.
We use a python module, cyvcf2, which quickly parses VCF files and puts information into objects that may be used to perform many operations beyond simply
calculating allele frequency differences between groups.
"""

import re
import sys
import os
from cyvcf2 import VCF
# see https://brentp.github.io/cyvcf2/docstrings.html

def filterVariantDepthMissData(gt_depths, gt_types, index2Group, Groups, Filter):
	PASS = 1
	sampleSizes_missing = {}
	for grp in Groups:
		# for each group, let's keep a list of the individuals that don't have sufficient sequencing depth,
		# as these individuals will be ignored when calculating allele frequencies.
		sampleSizes_missing[grp] = []
	for i in index2Group:
		if (gt_depths[i] < Filter['MinDepthPerInd']) or (gt_types[i] == 3): # has low depth or missing genotype
			grp = index2Group[i]
			sampleSizes_missing[grp].append(i)
	for grp in sampleSizes_missing:
		if len(sampleSizes_missing[grp]) > Filter['MaxMissIndPerCategory']:
			PASS = 0
	return(PASS, sampleSizes_missing)

def CalculateAFdiffByGroup(genotypes, index2Group, Groups, sampleSizes, sampleSizes_missing):
	AFbyGroup = {}
	GenosByGroup = {}
	for grp in Groups:
		AFbyGroup[grp] = 0
		GenosByGroup[grp] = [0,0,0] # Homozygous reference, heterozygous, homozygous alternate

	for i in index2Group:
		grp = index2Group[i]
		# need to ignore any missing individuals
		if i not in sampleSizes_missing[grp]:
			AFbyGroup[grp] += genotypes[i]
			GenosByGroup[grp][genotypes[i]] += 1
	for grp in Groups:
		# multiply denominator by 2 b/c you use genotypes to count AF, and these are diploid
		AFbyGroup[grp] = AFbyGroup[grp]/((sampleSizes[grp]-len(sampleSizes_missing[grp]))*2)
	return(AFbyGroup, GenosByGroup)

def main():
    if len(sys.argv) != 2:
        raise Exception('Must specify a VCF file')
    vcf_file_name = sys.argv[1]

    if not os.path.exists(vcf_file_name):
        raise Exception('The specified VCF file (%s) does not exist.' % vcf_file_name)

	vcf = VCF(vcf_file_name, gts012=True) # gts012=True makes value of 3 UNKNOWN, with 0,1,2 corresponding to numb ALT alleles

	'''
	Assign individuals to groups within which allele frequencies are computed
	I tried to code it such that groups may be removed, and those indices will not be iterated through
	'''
	Groups = []
	# 0: Tuskless MZ savanna elephants
	# 1: Tusked MZ savanna elephants
	# 2: Non-MZ savanna elephants (tusked)
	# 3: Forest elephants
	# 4: Asian elephants
	Groups.append(["0045B", "2983B", "2984B", "2986A", "G17A", "G19A", "G22A", "T2B"]) # T1A excluded, contaminated
	Groups.append(["2981B", "2982B", "2985B", "G18A", "G20A", "G21A"])
	Groups.append(["SRR958467", "SRR958468", "ERR2260496", "ERR2260497"])
	Groups.append(["Lcyc1", "Lcyc2"])
	Groups.append(["Emax1", "Emax2", "Emax3", "Emax4", "Emax5", "Emax6"])

	'''
	make list of samples to include or exclude in all downstream analyses
	here, we realized sample "T1A" was contaminated
	vcf.samples contains the samples we analyse downstream, and we can alter this list using the .set_samples method
	'''
	samplesToInclude = []
	samplesToExclude = ["T1A"]
	for sample in vcf.samples:
		if sample not in samplesToExclude:
			samplesToInclude.append(sample)
	vcf.set_samples(samplesToInclude) # only includes certain samples

	'''
	for each index in vcf.samples list that corresponds to an individual, which group does this individual belong to?
	while we're at it, let's collect the sample size per group, as this is used in the denominator to calculate allele frequency
	'''
	index2Group = {}
	sampleSizes = {}
	for grp in Groups:
		sampleSizes[grp] = 0
	for index in range(len(vcf.samples)):
		for grp in Groups:
			if vcf.samples[index] in Groups[grp]:
				index2Group[index] = grp
				sampleSizes[grp] += 1

	'''
	Specify minimum wsequencing deth and amount of missing data to tolerate
	'''
	FilterDP = {}
	FilterDP['MinDepthPerInd'] = 4 		# minimum sequencing depth per individual
	FilterDP['MaxMissIndPerCategory'] = 1	# minimum number of missing data/individuals per group

	# open output file and print the column names
	outFile = open('Results_allVariants.txt','w')
	col_names = [
		'Chrom',
		'Pos',
		'RefAllele',
		'AltAllele',
		'variantType',
		'Effect',
		'geneName',
		'geneNameEnsembl',
	]
	for grp in sorted(Groups):
		x = "AFGroup_" + str(grp)
		col_names.append(x)
	for grp in sorted(Groups):
		x = "GenosGroup_" + str(grp)
		col_names.append(x)
	print('\t'.join(col_names), file=outFile)

	# let's start parsing our VCF!
	for variant in vcf:
		# check if this site passed GATK VCF filtering, and only look at bi-allelic sites
		if variant.FILTER == "PASS" and len(variant.ALT) == 1:
			# check if this site has sufficient data, and keep track of how many individuals are missing per group as this number serves as denominator in calculating allele frequency
			PASS_depthMissData, sampleSizes_missing = filterVariantDepthMissData(variant.gt_depths, variant.gt_types, index2Group, Groups, FilterDP)

			if PASS_depthMissData == 1:
				print(variant.CHROM, " ", variant.end)
				# calculate allele frequency differences by group
				# NOTE: gt_types is array of 0,1,2,or 3, corresponding to HOM_REF, HET, HOM_ALT, and UNKNOWN, respectively
				AFbyGroup, GenosByGroup = CalculateAFdiffByGroup(variant.gt_types, index2Group, Groups, sampleSizes, sampleSizes_missing)
				# initialize variables... perhaps not necessary?
				effect = None
				geneName = None
				geneNameEnsembl = None
				for field in variant.INFO:
					if field[0] == 'ANN':
						info = field[1]
						# locus effect and gene name information separated by "|'
						infoList = info.split('|')
						effect = infoList[1]
						geneName = infoList[3]
						geneNameEnsembl = infoList[4]
				fields = [
					variant.CHROM,
					variant.end,
					variant.REF,
					variant.ALT,
					variant.var_type,
					effect,
					geneName,
					geneNameEnsembl,
				]
				for grp in sorted(Groups):
					fields.append(AFbyGroup[grp])
				for grp in sorted(Groups):
					fields.append(GenosByGroup[grp])
				print('\t'.join(fields), file=outFile)

if __name__ == '__main__':
    sys.exit(main())
