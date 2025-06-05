# Class 6: R functions
Leah Kim (A16973745)

- [1. Functions](#1-functions)
- [2. Generate DNA Function](#2-generate-dna-function)
- [3. Generate Protein Function](#3-generate-protein-function)

## 1. Functions

Let’s start writing our first silly function to add some numbers.

Every R function has 3 things - name (we get to pick this) - input
arguments (there can be loads of these separated by a comma) - the body
(the R code that does the work)

``` r
add <- function(x, y=10, z=0){
  x + y + z
}
```

I can just use this function

``` r
add(1, 100)
```

    [1] 101

``` r
add(1)
```

    [1] 11

Functions can either have “required” input arguments and “optional”
input arguments. The optional arguments are defined with an equals
default value (`y = 10`) in the function definition

``` r
add(1, 100, 10)
```

    [1] 111

## 2. Generate DNA Function

> Q. Write a function to return a DNA sequence of a user specified
> length. Call it `generate_dna()`

``` r
#generate_dna <- function(size = 5) {}

students <- c("jeff", "jeremy", "peter")
sample(students, size = 1)
```

    [1] "peter"

``` r
sample(students, size = 5, replace=TRUE)
```

    [1] "peter"  "jeremy" "peter"  "peter"  "jeremy"

Now work with `bases` rather than `students`

``` r
bases <- c("A", "C", "G", "T")
sample(bases, size = 10, replace = TRUE)
```

     [1] "A" "T" "A" "G" "A" "G" "C" "T" "A" "A"

Now I have a working ‘snippet’ of code I can use.

``` r
generate_dna <- function(size = 5) {
  bases <- c("A", "C", "G", "T")
sample(bases, size = size, replace = TRUE)
}
```

``` r
generate_dna(100)
```

      [1] "C" "G" "T" "A" "A" "C" "T" "T" "C" "C" "A" "G" "T" "T" "G" "T" "G" "A"
     [19] "A" "G" "C" "A" "A" "C" "C" "C" "G" "T" "A" "A" "G" "C" "G" "A" "T" "T"
     [37] "C" "A" "A" "G" "G" "C" "A" "T" "C" "G" "T" "C" "T" "T" "G" "C" "T" "G"
     [55] "C" "C" "G" "T" "A" "A" "T" "G" "G" "T" "A" "G" "T" "C" "G" "G" "G" "G"
     [73] "G" "T" "C" "T" "T" "C" "T" "T" "T" "C" "G" "C" "G" "T" "T" "T" "A" "A"
     [91] "G" "T" "G" "A" "A" "A" "C" "T" "A" "C"

``` r
generate_dna()
```

    [1] "T" "G" "C" "C" "C"

I want the ability to return a sequence like “AGTACCTG” i.e. a one
element vector where the bases are all together.

``` r
generate_dna <- function(size = 5, together = TRUE) {
  bases <- c("A", "C", "G", "T")
  sequence <- sample(bases, size = size, replace = TRUE)
  if(together) {
    sequence <- paste(sequence, collapse = "")
  } 
  return(sequence)
}
```

## 3. Generate Protein Function

We can ge thte set of 20 natural animo-acids from the **bio3d** package

> Q. Write a protein sequence generating function that will return
> sequences of a user specified length.

``` r
generate_protein <- function(size = 5, together = TRUE) {
  #get the 20 amino acids as a vector
  aa <- bio3d::aa.table$aa1[1:20]
  sequence <- sample(aa, size = size, replace = TRUE)
  
  #optionally return a single element string
  if(together) {
    sequence <- paste(sequence, collapse = "")
  } 
  return(sequence)
}
```

> Q. Generate random proteins equences of length 6-12 amino acids.

``` r
#generate_protein(6:12) returns an error
```

We can fix this inability to generate multiple sequences by either
editing and adding to the function body code (e.g. a for loop) or by
using the R **apply** family of utility functions

``` r
ans <- sapply(6:12, generate_protein)
ans
```

    [1] "FRQGKV"       "WDLTQYA"      "RQHAEVQP"     "KLIDFWDEL"    "VPMGFQAFIQ"  
    [6] "MFPMCALPYLD"  "WRMPVGDCAICI"

It would be cool and useful if I could get FASTA format output. I want
this to look like

    >ID.6
    HLDVLV
    >ID.7
    VREAIQN
    >ID.8
    WPRSKACN

The functions `cat` and `paste` can help us here.

``` r
ans <- sapply(6:12, generate_protein)
cat(paste(">ID.",6:12, sep="", "\n", ans), sep ="\n")
```

    >ID.6
    WPHYWP
    >ID.7
    HLGIVMQ
    >ID.8
    NLVCEKNG
    >ID.9
    SPPRDKDDH
    >ID.10
    MTKPRTAWIR
    >ID.11
    EWENRGCSSNK
    >ID.12
    SSQWNECMFDPD

> Q. Determine if any of these sequences can be found in nature or are
> they unique? Why or why not?

The sequences generated are as followes:

    >ID.6
    SFDRHS
    >ID.7
    HQNLFYY
    >ID.8
    DSMEMNDL
    >ID.9
    STYFCEKGC
    >ID.10
    CVDIIEFNKR
    >ID.11
    ECFMCPHRVDN
    >ID.12
    PSRKPESIFEHE 

I BlastP my FASTA format sequences against NR and found that the
sequences of lengths 6 through 8 are not unique and found in the
databases with 100% coverage and identity.

Random sequences of length 9 and above are unique and cannot be foud in
the databases.
