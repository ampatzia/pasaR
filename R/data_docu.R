#Datasets
# 1) Streptococcus and Bacillus pangenomic data
# 2) Streptococcus and Bacillus pangenomic data names


#' Streptococcus and Bacillus pangenomic data
#'
#' A cluster membership data frame for genes from 81 bacterial genomes (rows) spread
#' in 149721 clusters (columns).
#' The dataset was created from publicly available data from Ensembl,It was produced
#' by a standard clustering pipeline with the default settings (BLAST and MCL).
#' and contains seven (81) strains of the following bacteria :
#' twelve (12) of Streptococcus pneumoniae, thirteen (13) of Streptococcus
#' Pyogenes, thirty nine (39) of Bacillus cereus and seventeen (17) of
#' Bacillus thuringiensis. Genome names are contained in the bac_stre_names
#' dataset.
#'
#' @usage data(bac_stre)
#'
#'
#' @references Mpatziakas A, Psomopoulos FE, Moysiadis T and Sgardelis S. Computing pangenome statistics in R
#'  F1000Research 2017, 6(ISCB Comm J):1529 (poster) (doi: 10.7490/f1000research.1114765.1)
#'
#'
#' @examples
#' data(bac_stre)
"bac_stre"


#' Names of Streptococcus and Bacillus species
#'
#' This dataset contains the names of the genomes inside the bac_stre dataset.
#'
#'
#' @usage data(bac_stre_names)
#'

#'
#' @references Mpatziakas A, Psomopoulos FE, Moysiadis T and Sgardelis S. Computing pangenome statistics in R
#'  F1000Research 2017, 6(ISCB Comm J):1529 (poster) (doi: 10.7490/f1000research.1114765.1)
#'
#'
#' @examples
#' data(bac_stre_names)
"bac_stre_names"
