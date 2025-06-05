# Class 17: Analyzing sequencing data in the cloud
Leah Kim A16973745

- [Setup](#setup)
- [PCA](#pca)
  - [Plot of Transcriptomic Profiles](#plot-of-transcriptomic-profiles)

## Setup

``` r
library(tximport)

# setup the folder and filenames to read
folders <- dir(pattern="SRR21568*")
samples <- sub("_quant", "", folders)
files <- file.path( folders, "abundance.h5" )
names(files) <- samples

txi.kallisto <- tximport(files, type = "kallisto", txOut = TRUE)
```

    1 2 3 4 

``` r
head(txi.kallisto$counts)
```

                    SRR2156848 SRR2156849 SRR2156850 SRR2156851
    ENST00000539570          0          0    0.00000          0
    ENST00000576455          0          0    2.62037          0
    ENST00000510508          0          0    0.00000          0
    ENST00000474471          0          1    1.00000          0
    ENST00000381700          0          0    0.00000          0
    ENST00000445946          0          0    0.00000          0

Number of transcripts for each sample:

``` r
colSums(txi.kallisto$counts)
```

    SRR2156848 SRR2156849 SRR2156850 SRR2156851 
       2563611    2600800    2372309    2111474 

How many transcripts detected in at least one sample:

``` r
sum(rowSums(txi.kallisto$counts)>0)
```

    [1] 94561

Filter out transcripts with no reads:

``` r
to.keep <- rowSums(txi.kallisto$counts) > 0
kset.nonzero <- txi.kallisto$counts[to.keep,]
```

Filter out transcripts with no change:

``` r
keep2 <- apply(kset.nonzero,1,sd)>0
x <- kset.nonzero[keep2,]
```

## PCA

Run PCA

``` r
pca <- prcomp(t(x), scale=TRUE)
summary(pca)
```

    Importance of components:
                                PC1      PC2      PC3   PC4
    Standard deviation     183.6379 177.3605 171.3020 1e+00
    Proportion of Variance   0.3568   0.3328   0.3104 1e-05
    Cumulative Proportion    0.3568   0.6895   1.0000 1e+00

### Plot of Transcriptomic Profiles

PC1 vs PC2:

``` r
plot(pca$x[,1], pca$x[,2],
     col=c("blue","blue","red","red"),
     xlab="PC1", ylab="PC2", pch=16)
```

![](class17_files/figure-commonmark/unnamed-chunk-8-1.png)

PC1 vs PC3:

``` r
plot(pca$x[,1], pca$x[,3],
     col=c("blue","blue","green","green"),
     xlab="PC1", ylab="PC3", pch=16)
```

![](class17_files/figure-commonmark/unnamed-chunk-9-1.png)

PC2 vs PC3:

``` r
plot(pca$x[,2], pca$x[,3],
     col=c("red","red","green","green"),
     xlab="PC2", ylab="PC3", pch=16)
```

![](class17_files/figure-commonmark/unnamed-chunk-10-1.png)
