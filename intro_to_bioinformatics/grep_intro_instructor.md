# Genome tips and tricks

## Manipulating files with grep
### What is grep?

grep is a powerful command-line search tools that is included as part of most Unix-like systems. It is one of the most useful tools in bioinformatics! At the most basic level, grep searches for a string of characters that match a pattern and will print lines containing a match. Basic syntax is:  

`grep 'pattern' file_to_search`  


By default, grep will match any part of the string, so for example:

```
echo 'my dog is brown' > sample.txt  
grep 'dog' sample.txt  
grep 'do' sample.txt  
grep 'd' sample.txt
```


The above three grep commands will all match the text in the example file and will print the line. However, there are a huge number of arguments that can modify how grep behaves. Here are a few useful examples!


`grep -w` matches *entire words*.
- So in the above example:  
`grep -w 'dog' sample.txt` would match the string, but `grep -w 'do' sample.txt` would not.

`grep -i` allows case-insensitive matches
- In the above example, `grep -i 'DOG'` would still match the line

`grep -v` *inverts*, returning lines that *do not* match the pattern.

`grep -o` returns only the matching words, not the entire line.

`grep -c` counts the number of lines that match the pattern.
- Equivalent to `grep 'pattern' file | wc -l`

Print lines before/after a match:
- `grep -A [n]` returns matching line and *n* lines after match
- `grep -B [n]` returns matching line and *n* lines before match
- `grep -C [n]` returns matching line and *n* lines before *and* after match

Note: `grep -A` is very useful for pulling out specific lines from a FASTA file...just make sure your FASTA file is single line and not multi-line!

`grep -f pattern.txt` matches a list of patterns contained in pattern.txt against target file.
- I.e. the command `grep -f patterns.txt file.txt` will match all patterns in patterns.txt against file.txt
- Pattern file has one pattern per line


Remember that we can combine different grep arguments with each other! E.g.:

`grep -w -A 1 '>X' data/dmel-all-chromosome-r6.20.fasta`  
What will this line do??

`grep -w -v -f first3.txt data/human_enhancers.bed | less`  
How about this line?

There are many other functions of grep! When in doubt, remember you can check the help page using `man grep`...or just by using google :)

### Pattern matching with regular expressions
Regular expressions, aka "regex", are patterns that describe sets of strings. In other words, they allow you to match complex patterns with grep, not just exact matches. Regex is extremely powerful, but can also get (very) complicated, so we'll just stick to a few basic uses.

Regex has certain meta-characters that are reserved for special uses:
```
^: matches pattern at start of string
$: matches pattern at end of string
.: matches any character except new lines
[]: matches any of enclose characters
[^]: matches any characters *except* ones enclosed (note: is different from ^)
\: "escapes" meta-characters, allows literal matching
```

Note that if we want to match any of these special characters literally (e.g. matching a period character "."), we would need to use a "\\" to escape it first:  
`grep '\.' data/dmel-all-no-analysis-r6.20.gff` will literally match a period.

One example of how regex can come in handy is using the ^ special character to quickly count how many sequences are in a FASTA file, which we would do as follows:

`grep -c '^>' data/mel-all-chromosome-r6.20.fasta`

This command matches lines in the FASTA file that start with a ">" character, i.e. the header lines, and uses the -c argument to count how many matches!


### Practice:

- How many sequences are in the file Falb_COL2.1.fastq?

`grep -c '^@' data/Falb_COL2.1.fastq`

- Filter out lines matching unassigned contigs (chrUn) in the file hg19.genome and direct the output to a file.

`grep -v 'chrUn' data/hg19.genome  > hg19_noUn.genome`

- Use grep to pull out the header line of **only** the 2R chromosome arm from dmel-all-chromosome-r6.20.fasta and direct the output to a file.

`grep -w '^>2R' data/dmel-all-chromosome-r6.20.fasta > 2R_header.txt`

- Use grep to extract the lines of only the *major chromosome arms* (2L, 2R, 3L, 3R, and X) from the file dmel-all-no-analysis-r6.20.gff and pipe to less

```
printf '^>%s\n' 2L 2R 3L 3R X > major_arms.txt
grep -w -f major_arms.txt \
data/dmel-all-no-analysis-r6.20.gff | less
```
