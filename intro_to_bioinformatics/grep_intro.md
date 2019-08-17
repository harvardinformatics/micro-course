# Genome tips and tricks

## Manipulating files with grep
### What is grep?

grep is a powerful command-line search tools that is included as part of most Unix-like systems. It is one of the most useful tools in bioinformatics! At the most basic level, grep searches for a string of characters that match a pattern and will print lines containing a match. Basic syntax is:  

`grep "pattern" file_to_search`  


By default, grep will match any part of the string, so for example:

```
echo "my dog is brown" > sample.txt  
grep "dog" sample.txt  
grep "do" sample.txt  
grep "d" sample.txt
```


The above three grep commands will all match the text in the example file and will print the line. However, there are a huge number of arguments that can modify how grep behaves. Here are a few useful examples!

- `grep -v` *inverts*, returning lines that *do not* match the pattern.

`grep -w` matches *entire words*.
- So in the above example:  
`grep -w "dog" sample.txt` would match the string, but `grep -w "do" sample.txt` would not.

`grep -o` returns only the matching words, not the entire line.

`grep -c` counts the number of lines that match the pattern.
- Equivalent to `grep "pattern" file | wc -l`

Print lines before/after a match:
- `grep -A [n]` returns matching line and *n* lines after match
- `grep -B [n]` returns matching line and *n* lines before match
- `grep -C [n]` returns matching line and *n* lines before *and* after match

Note: `grep -A` is very useful for pulling out specific lines from a FASTA file...just make sure your FASTA file is single line and not multi-line!

`grep -f pattern_file.txt` matches a list of patterns contained in pattern_file.txt against target file.


### Practice:

- How many sequences are in the file Falb_COL2.1.fastq?

`grep -c "^@" Falb_COL2.1.fastq`

- Filter out lines matching unassigned contigs (chrUn) in the file hg19.genome

`grep -v "chrUn" hg19.genome | less`

- Use grep to pull out the header line of **only** the 2R chromosome arm from dmel-all-chromosome-r6.20.fasta

- Use grep to extract the lines of only the *major chromosome arms* (2L, 2R, 3L, 3R, and X) from the file dmel-all-no-analysis-r6.20.gff
