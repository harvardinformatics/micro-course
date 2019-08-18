# Genome tips and tricks

## Manipulating files with awk
### What is awk?

Invented in the 1970's, awk is a scripting language included in most Unix-like operating systems. It specializes in one-liner programs and manipulating text files.

In many cases, if you're parsing information from a text file (such as a BED file, FASTA file, etc.), you could write a Python script...or you could do it with awk in a single line!

### Syntax
awk scripts are organized as:

`awk 'pattern { action; other action }' file`

Meaning that every time that the pattern is true, awk will execute the action in the brackets. By default the pattern matches every line, so the action will be taken every line in the input file, e.g. the following command that prints every line:

`awk '{print}' data/hg38.genome | less`

The two most important patterns are `BEGIN` and `END`, which tell the action to take place before any lines are read and after the last line.

 `awk 'BEGIN{sum=0} {sum+=1} END {print sum}' data/hg38.genome`

 The above line sets a variable at the start of the script, adds 1 to it every line, then prints its value at the end.


### Input and output
Input to awk is split into **records** and **fields**.
- By default, **records** are separated by newline character, i.e # of records = # of lines in input file
- Each record is subdivided into **fields**, i.e. columns, as determined by the field separator (see below)

There are several important built-in variable in awk. The fields (columns) of each record are referred to by `$number`, so the first column would be `$1`, second would be `$2`, etc. `$0` refers to the entire record.<br/>
So to print the second column of each line in the file, we'd use:

`awk '{print $2}' data/hg38.genome | less`

And if we wanted to print the second then the first:

`awk '{print $2,$1}' data/hg38.genome | less`

Note that the different fields are joined with commas when printing.
<br/>
<br/>

We can also print strings using using quotation marks:

`awk '{print "First column:" $1}' data/hg38.genome | less`

Which for every line of the file will print the text "First column:" followed by the value in the first field.

---
awk has several other built-in variables that are very useful for parsing text:

>FS: field separator (default: white space)<br/>
OFS: output field separator, i.e. what character separates fields when printing<br/>
RS: record separator, i.e. what character records are split on (default: new line)<br/>
ORS: output record separator<br/>
NR: number of records in input (# lines by default)

Using these, we can convert between file formats, e.g. make a comma-separated text file into a tab-separated file:

`awk 'BEGIN{FS="," ; OFS="\t"} {print $0}' data/enhancers.csv > data/enhancers.tsv`


### Conditionals and pattern matching
Like other programming languages, awk allows conditional matching with if/else statements.

`awk '{if(condition) action; else other_action}'`

awk uses the following conditional operators:

>'==': equal to  
'!=': not equal to  
'>': greater than  
'>=': greater than or equal to  
'<': less than    
'<=': less than or equal to  
'&&': AND  
'||': OR

In addition, awk also supports string matching using regular expressions, using the following expressions:

>'\~': matches  
'!~' does not match

For string matching, the pattern being matched must be enclosed by slashes, like so:

`awk '{if($1 ~ /pattern/) print}'`

---
### Example uses:
- Only print annotations on a specific scaffold (chr1) that fall between 1Mb and 2Mb:

`awk 'BEGIN{FS="\t";OFS="\t"} {if($1 == "chr1" && $2 >=1000000 && $2 <= 2000000) print}' data/ucscGenes.bed | less`

**Note**: when we specify that we only want annotations from chr1, we're using exact match (`== "chr1"`) and not pattern match (`~ /chr1/`)...why is this??

- Only print lines of GFF file that match the string "exon" in their third column:

`awk 'BEGIN{FS="\t"} {if($3 ~ /exon/) print $0}' data/dmel-all-no-analysis-r6.20.gff | less`

- Convert from GFF (genome feature file) to BED file

`grep -v '^#' data/dmel-all-no-analysis-r6.20.gff | awk 'BEGIN{FS="\t"; OFS="\t"} {print $1,$4-1,$5}' | less`

**Note**: remember that BED and GFF files have different coordinate systems, i.e. BED start coordinate is 0 based, half-open, GFF is 1-based inclusive! Also, we are first using grep to skip the header lines in the GFF file.


### Practice
Using awk:<br/>
- Pull out only the CDS annotations from the GFF file dmel-all-no-analysis-r6.20.gff and output them in BED format

`awk 'BEGIN{FS="\t"; OFS="\t"} {if($3 ~ /CDS/) print $1,$4-1,$5}' data/dmel-all-no-analysis-r6.20.gff`

- Extract FASTA information from the SAM file Falb_COL2.final.sam (hint: you will need to remove the header lines first, they start with @!) (Another hint: sequence ID = 1st column, sequence = 10th column)

`grep -v '^@' data/Falb_COL2.final.sam | awk 'BEGIN{FS="\t"; OFS="\n"} {print ">"$1,$10}' | less`

- Write a command to calculate that average of the 5th column (i.e. mapping quality score) of a tab-separated SAM file Falb_COL2.final.sam and output it (hint: as above, need to remove headers)

`grep -v '^@' data/Falb_COL2.final.sam | awk 'BEGIN{FS="\t"; sum=0} {sum+=$5} END{print sum/NR}'`

- Calculate the average length of gene annotations *only on the 2L arm* from the file dmel-genes.bed (hint: you'll need to use a combination of grep and awk...)

`grep "2L" data/dmel-genes.bed | awk 'BEGIN{FS="\t"; sum=0} {len=$3-$2; sum=sum+len} END{print sum/NR}'`
