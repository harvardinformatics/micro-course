# A small bioinformatics exercise
This notebook uses a small bioinformatics exercise to show aspects of the Python programming 
language in the context of a real(ish) data processing activity.

We will be reading, writing, and manipulating text files and running a small sequence alignment
program.  Over the course of this we will cover programming topics such as:

   * Built-in Python types including strings, ints, floats
   * Python code blocks including if/then/else, while and for loops, functions and classes,
     and context managers
   * Exception handling
   * Context managers
   * System calls
   
## Five ways to setup a filename
This section shows five different ways to get to a filename that can be opened.  At the same time
we will go over a number of basics about Python including it's object-oriented types,
string formatting and a little bit about lists.

### 1. Assign a string literal to a variable
You can use tab completion to fill out the filename.  Python will automatically assign a
type to the variable

### 2. Concatenate string elements.  Try it as a function
Strings can be concatenated with the '+' operator.  Non-strings must be
converted first.

#### Define a function with positional and keyword arguments.  Don't forget to indent.

#### Incorrect indentation is a syntax error in Python

#### Call the function with just positional arguments

#### Call the function with reordered arguments as keywords

### 3. Formatted strings
Python supports both positional and named string template substitution.  See the
[Pyformat page](https://pyformat.info/) for details

#### String concatentation is expensive because Python strings are immutable
Because strings cannot be changed, string concatenation requires a new string
to be created for each concatenation operation.

#### Old style string formatting using formatting codes is common
Old style formatting uses format codes like _%s_ and _%d_ as place holders in
a string template.

#### format function is more readable and powerful
The format function of strings allows for positional substitution like old style
formatting, but also supports named place holders and rich formatting options

### A brief interlude about classes, functions, and objects in Python



### 4. Joining list elements
A list of elements can be joined.  This is particularly useful when you have a loop that 
accumulates a variable number of values.

#### Like arrays in other languages, Python lists are a group of items that can be indexed by an integer

### 5. Joining list elements with os.path.join
The _os_ module must be imported and contains functions that are sensitive to the operating system

## Interacting with files
This section shows how to read, write and process files

