# Genome tips and tricks

## Manipulating files with grep
### What is grep?

grep is a powerful command-line search tools that is included as part of most Unix-like systems. It is one of the most useful tools in bioinformatics! At the most basic level, grep searches for a string of characters that match a pattern and will print lines containing a match. Basic syntax is:  

`grep "pattern" file_to_search`  


By default, grep will match any part of the string, so for example:

`echo "my dog is brown" > sample.txt`  
`grep "dog" sample.txt`  
`grep "do" sample.txt`  
`grep "d" sample.txt`  

The above three grep commands will all match the text in the example file and will print the line. However, there are a huge number of arguments that can modify how grep behaves
