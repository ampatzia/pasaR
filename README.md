[![Build Status](https://travis-ci.org/ampatzia/pasaR.svg?branch=master)](https://travis-ci.org/ampatzia/pasaR) [![Rhub Status](https://img.shields.io/badge/R--hub%20builder-Ok-brightgreen.svg)](https://builder.r-hub.io/status/original/pasaR_0.9.0.tar.gz-15185b02e6b84de9b3103203db732b84)
[![Coverage Status](https://img.shields.io/codecov/c/github/ampatzia/pasaR/master.svg)](https://codecov.io/github/ampatzia/pasaR?branch=master)

# R package pasaR
We present an R package, named pasaR, usable in the later stages of an pangenomic
analysis, i.e. after the construction of the gene families for a given set of genomes, based on information of the full complement of gene families. A complete methodology is proposed, suitable for sets of genomes of varying complexity, optimizing and enriching an assortment of existing measures from micropan, the only R package currently available on CRAN for such studies. This is an on-going project so better documentantion, a more extensive vignette and additional functions will be added. However the package is fully functional.

# Install package

If package *devtools* is present, simply run:

```R
library(devtools)
install_github("ampatzia/PasaR")
```

# Vignettes

Currently there are available two vignettes for the package, both precompiled in pdf:

* [Small case study and example of usage in a pangenome of 81 bacterial strains.](https://github.com/ampatzia/pasaR/blob/master/vignettes/Pangenome_analysis_with_pasaR.pdf) .
* [Benchmark comparison to package micropan](https://github.com/ampatzia/pasaR/blob/master/vignettes/Benchmark_Comparison.pdf) .

# Relevant publications

Mpatziakas A, Psomopoulos FE, Moysiadis T and Sgardelis S. Computing pangenome statistics in R. F1000Research 2017, 6(ISCB Comm J):1529 (poster) ([doi: 10.7490/f1000research.1114765.1](http://dx.doi.org/10.7490/f1000research.1114765.1))

# Contact
Any questions should be directed to ampatziakas at gmail.com
