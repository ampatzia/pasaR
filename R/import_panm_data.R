#' @keywords  internal
globalVariables(".")


#'  Make panmatrix (MCL data)
#'
#' This function allows importing MCL output
#' @param file_path Disk Path to file
#' @export
#' @examples \dontrun{make_panmatrix(file_path)}
#' @note  MCL output as described in F. E. Psomopoulos, O. T. Vrousgou, and P. A. Mitkas, "Large-scale modular comparative genomics: the Grid approach [v1; not peer reviewed]," F1000research 2015, vol. 4(ISCB Com, iss. 377, p. 1, 2015. doi:10.7490/f1000research.1110127.1
#'        A. M. Kintsakis, F. E. Psomopoulos, and P. A. Mitkas, "Data-aware optimization of bioinformatics workflows in hybrid clouds," Journal of big data, vol. 3, iss. 20, pp. 1-26, 2016. doi:10.3389/fpls.2016.00554
#' @importFrom magrittr "%>%"
#' @importFrom ggplot2 geom_bar geom_line geom_point ggtitle scale_y_continuous xlab ylab aes


make_panmatrix <- function(file_path){

  Organism<-NULL
  Protein<-NULL
  Cluster<-NULL
  Proteins<-NULL
  cluster_composition<-NULL
  x1<-NULL

  make_base_df <- function(x){
    work_list <- scan(file=x,what="character,n=195,", sep=" ", allowEscapes = TRUE)%>%
      stringr::str_split_fixed( ., " ", n = Inf) %>%
      sapply(., stringi::stri_escape_unicode) %>%  #Escapes all Unicode (not ASCII-printable) code points ie. single /
      sapply(., function(x) stringr::str_split_fixed(x, "[\\\\]+t|[^[:print:]]" , n = Inf)) %>%
      lapply(., function(x) stringr::str_split_fixed(x, "\\|", n=Inf))%>%
      lapply(., function(x) data.frame(x, stringsAsFactors=FALSE)) %>%
      lapply(., function(x){colnames(x)[1] <- "x1"; x}) %>%
      lapply(., function(x) tidyr::separate(x,x1,into = c("Organism", "Protein", "Other"), sep="\\$"))%>%
      lapply(., function(x) x[,names(x) %in% c("Organism", "Protein") ])

    for (i in 1:length(work_list)){

      work_list[[i]]<-transform(work_list[[i]],Cluster=i)
    }


    #Make list to dataframe
    result_df <- dplyr::bind_rows(work_list)
    rm(work_list)

    return(result_df)

  }

  cluster_composition <- function(x){
    result_df_cl1<-x %>% dplyr::group_by(., Cluster,Organism) %>%
      dplyr::summarise(., Proteins=length(Protein))
    return(result_df_cl1)}


  panm <- file_path %>% make_base_df(.) %>% cluster_composition(.) %>%tidyr::spread(.,Cluster,Proteins,fill=0)
  #organism_names<-panm[,1]
  panm<-panm[,-1]
  return(panm)

}

#' Make panmatrix (fami MCL data)
#'
#' This function allows importing MCL output
#' @param file_path Disk path to file
#' @export
#' @examples \dontrun{make_panmatrix_fami(file_path)}
#'

make_panmatrix_fami<-function(file_path){
  Organism<-NULL
  Protein<-NULL
  Cluster<-NULL
  Proteins<-NULL
  cluster_composition<-NULL
  V1<-NULL

  work_list <- readr::read_delim(file_path, "\t", escape_double = FALSE, col_names = FALSE,
                                 trim_ws = TRUE)

  split <- as.data.frame(stringr::str_split_fixed(work_list$X1," ", n=Inf))

  work_list <- dplyr::bind_cols(split,work_list[,-1])
  work_list$V2 <- as.character(work_list$V2)
  work_list <- tidyr::gather(work_list,V1)
  work_list <- work_list[, c(1,3)]
  colnames(work_list) <- c("Cluster", "V1")
  work_list <- work_list[stats::complete.cases(work_list),]
  work_list <- tidyr::separate(work_list, V1, into=c("n1", "n2", "version", "PID"), sep="-")
  work_list$n1 <- paste0(work_list$n1, "_", work_list$n2)
  work_list$Cluster <- as.numeric(unlist(stringr::str_extract_all(work_list$Cluster, "\\(?[0-9,.]+\\)?")))

  work_list <- work_list[, c(1,2,5)]
  colnames(work_list)<-c("Cluster", "Organism", "Protein")

  cluster_composition <- function(x){
    result_df_cl1 <- x %>% dplyr::group_by(., Cluster, Organism) %>%
      dplyr::summarise(., Proteins=length(Protein))
    return(result_df_cl1)
  }



  panm <- work_list %>% cluster_composition(.) %>%tidyr::spread(.,Cluster,Proteins,fill=0)

  org_names <- panm[,1]
  panm <- panm[,-1]
  }


#' Make panmatrix (fami 2 MCL data)
#'
#' This function allows importing BLAST MCL output with default parameters
#' @param file_path file Path to file
#' @export
#' @examples \dontrun{make_panmatrix_fami2(file_path)}
#'

make_panmatrix_fami2 <- function (file_path){
  Organism<-NULL
  Protein<-NULL
  Cluster<-NULL
  cluster<-NULL
  Proteins<-NULL
  cluster_composition<-NULL
  value<-NULL
  V1<-NULL

  work_list <- readr::read_delim(file_path, "\t", escape_double = FALSE,
                                 col_names = FALSE, trim_ws = TRUE)

  work_list$V1 <- (paste0("cluster", 1:nrow(work_list)))
  colvals <- paste0("X", 1:(ncol(work_list)-1))
  work_list <- tidyr::gather(work_list, cluster,value, 1:(ncol(work_list)-1))
  work_list <- work_list[, c(1, 3)]
  colnames(work_list) <- c("Cluster", "V1")
  work_list <- work_list[stats::complete.cases(work_list), ]
  work_list <- tidyr::separate(work_list, V1, into = c("n1", "n2",
                                                       "version", "PID"), sep = "-")
  work_list$n1 <- paste0(work_list$n1, "_", work_list$n2)
  work_list$Cluster <- as.numeric(unlist(stringr::str_extract_all(work_list$Cluster,
                                                                  "\\(?[0-9,.]+\\)?")))
  work_list <- work_list[, c(1, 2, 5)]
  colnames(work_list) <- c("Cluster", "Organism", "Protein")
  cluster_composition <- function(x) {
    result_df_cl1 <- x %>% dplyr::group_by(., Cluster, Organism) %>%
      dplyr::summarise(., Proteins = length(Protein))
    return(result_df_cl1)
  }
  panm <- work_list %>% cluster_composition(.) %>% tidyr::spread(.,
                                                                 Cluster, Proteins, fill = 0)
  org_names <- panm[, 1]
  panm <- panm[, -1]
}
