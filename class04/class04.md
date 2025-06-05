# Intro to R
Leah Kim A16973745
2025-04-10

``` r
#My first R script
x <- 1:50

plot(x, sin(x))
```

![](class04_files/figure-commonmark/unnamed-chunk-1-1.png)

``` r
plot(x, sin(x), typ="l", col="red",lwd=3, xlab="Silly x axis", ylab="Sensible y axis")
```

![](class04_files/figure-commonmark/unnamed-chunk-1-2.png)
