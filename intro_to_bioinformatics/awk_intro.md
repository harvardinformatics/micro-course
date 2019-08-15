# Genome tips and tricks

## Manipulating files with awk
### What is awk?

Invented in the 1970's, awk is a scripting language included in most Unix-like operating systems. It specializes in one-liner programs and manipulating text files.

In many cases, if you're parsing information from a text file (such as a BED file, FASTA file, etc.), you could write a Python script...or you could do it with awk in a single line!

### Syntax
awk scripts are organized as:

`pattern { action; other action }`

Meaning that every time that the pattern is true, awk will execute the action in the brackets. By default the pattern matches every line, so the action will be taken every line, e.g. the following command that prints every line:

`awk '{print}' input.txt`

The two most important patterns are `BEGIN` and `END`, which tell the action to take place before any lines are read and after the last line.

 `awk 'BEGIN{sum=0} {sum+=1} END{print sum}'`

 The above line sets a variable at the start of the script, adds 1 to it every line, then prints its value at the end.


### Input and output
Input to awk is split into **records** and **fields**.
- By default, **records** are separated by newline character, i.e # of records = # of lines in input file
- Each record is subdivided into **fields**, i.e. columns, as determined by the field separator (see below)

There are several important built-in variable in awk. The fields (columns) of each record are referred to by `$number`, so the first column would be `$1`, second would be `$2`, etc. `$0` refers to the entire record.<br/>
So to print the fifth column of each line in the file, we'd use:

`awk '{print $5}'`

And if we wanted to print the first, third, and sixth:

`awk '{print $1,$3,$6}'`

Note that the different fields are joined with commas when printing.
<br/>
<br/>

We can also print strings using using quotation marks:

`awk '{print "First column:" $1}'`

Which for every line of the file will print the text "First column:" followed by the value in the first field.

---
awk has several other built-in variables that are very useful for parsing text:

>FS: field separator (default: white space)<br/>
OFS: output field separator, i.e. what character separates fields when printing<br/>
RS: record separator, i.e. what character records are split on (default: new line)<br/>
ORS: output record separator<br/>
NR: number of records in input (# lines by default)

Using these, we can convert between file formats, e.g. make a comma-separated text file into a tab-separated file:

`awk 'BEGIN{FS="," ; OFS="\t"} {print $0}'`


### Conditionals and pattern matching
Like other programming languages, awk allows conditional matching with if/else statements.

`awk '{if(condition) action; else other_action'}`

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
Here are some example uses:
- Only display hits in a tab-separated BLAST file with at least 90% identity to the query sequence:

`awk 'BEGIN{FS="\t"} {if($4 >= 90) print $0}' file.blast`

- Only print lines of VCF file that match the string "chr" in their first column:

`awk 'BEGIN{FS="\t"} {if($1 ~ /chr/) print $0}' file.vcf`

- Convert from GFF (genome feature file) to BED file

`awk 'BEGIN{FS="\t"; OFS="\t"} {print $1,$3,$4}' file.gff`

### Practice
Using awk:<br/>
- Subset a BED file to only include lines between certain coordinates
- Extract FASTA information from a SAM file
- Write a command to calculate that average of the third column of a tab-separated list and output it.

>Hints:
- FASTA format: <br/>
\>ID line<br/>
ATGCGCGTCAA...<br/>
