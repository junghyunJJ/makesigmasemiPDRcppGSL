# makesigmasemiPDRcppGSL
The makesigmasemiPDRcppGSL was r package for [marsR](https://github.com/junghyunJJ/marsR) to compute the positive semidefinite of LD matrix.


## Installation
- The C++ library for GNU [GSL](https://www.gnu.org/software/gsl/) is required.
- makesigmasemiPDRcppGSL works only on ***nix (Linux, Unix such as macOS) system**. please check **.Platform$OS.type** function.
- We currently only support R 3.5+.*

```
install.packages("Rcpp")
install.packages("RcppGSL")
install.packages("devtools")
devtools::install_github("junghyunJJ/makesigmasemiPDRcppGSL")
```

## Example
```
> library(makesigmasemiPDRcppGSL)
> data("testdata")
> dim(testdata$geno)
[1] 471 489
> testdata$geno[1:5,1:5]
    sample1 sample2 sample3 sample4 sample5
rs1       1       1       1       1       2
rs2       0       0       0       0       0
rs3       0       1       0       0       0
rs4       0       1       0       0       0
rs5       0       1       0       0       0
> res <- makesigmasemiPDRcppGSL(testdata$geno)
> dim(res)
[1] 471 471
> res[1:5,1:5]
            [,1]        [,2]        [,3]        [,4]        [,5]
[1,]  1.10000000  0.31751143 -0.30376036  0.09944961  0.09944961
[2,]  0.31751143  1.10000000 -0.07018644 -0.08088798 -0.08088798
[3,] -0.30376036 -0.07018644  1.10000000 -0.03654935 -0.03654935
[4,]  0.09944961 -0.08088798 -0.03654935  1.10000000  1.00000000
[5,]  0.09944961 -0.08088798 -0.03654935  1.00000000  1.10000000
```

