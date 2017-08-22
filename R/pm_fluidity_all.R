#'  Genome Fluidity
#'
#' This function computes fluidity without sampling.
#'
#'
#' @param Panmatrix Panmatrix produced by make_panmatrix functions
#' @export
#' @details  Fluidity takes values in [0,1], with 1 denoting no common genes.This metric was introduced is called the Sorensen distance (Deza & Deza, 2009) and was first introduced in the
#'   context of a pangenome analysis in (Kislyuk et al ,2011).
#'
#' @examples \dontrun{pm_fluidity_all(Panmatrix)}
#'
#' @references
#' A. O. Kislyuk, B. Haegeman, N. H. Bergman, and J. S. Weitz, "Genomic fluidity:
#'  an integrative view of gene diversity within microbial populations," BMC genomics, pp. 12-32, 2011.
#' M. M. Deza and E. Deza, Encyclopedia of Distances. Springer, 2009.

pm_fluidity_all <- function (Panmatrix){
  Genome_1<-NULL #prevent namespace problems
  Genome_2<-NULL
  Fluidity<-NULL

  all_comb <- as.data.frame(gtools_comb(nrow(Panmatrix),2))
  panm <- sapply(Panmatrix, function(x) as.logical(x))

  fluid <- function(x){
    g1 <- panm[x,] %>% .[,colSums(.)>0]
    flu <- sum(abs(g1[1,]-g1[2,]))/sum(colSums(g1))

    return(flu)}

  all_comb$fluidity <- apply(all_comb,1,fluid)
  colnames(all_comb) <- c("Genome_1","Genome_2","Fluidity")



  test1 <- split(all_comb,all_comb$Genome_1)
  test2 <- split(all_comb,all_comb$Genome_2)
  dummy <- data.frame(Genome_1=as.numeric(NA),Genome_2=(NA),Fluidity=(NA))
  test1[[max(all_comb$Genome_2)]] <- dummy
  names(test1)[max(all_comb$Genome_2)] <- as.character(max(all_comb$Genome_2))


  test2$'1' <- dummy
  test1 <- test1[order(as.numeric(names(test1)))]
  test2 <- test2[order(as.numeric(names(test2)))]


  test2 <- lapply(test2, function(x){x<-x[,c(2,1,3)]})
  test2 <- lapply(test2,stats::setNames,colnames(test1$`1`))

  tester_merge <- Map(rbind,test2, test1)


  all_fluid <- do.call(rbind,tester_merge)%>%dplyr::arrange(.,Genome_1,Genome_2)
  all_fluid <- all_fluid[stats::complete.cases(all_fluid),]
  res <- all_fluid%>% dplyr::group_by(.,Genome_1)%>%dplyr::summarise(.,Fluidity=mean(Fluidity))%>%dplyr::arrange(.,dplyr::desc(Fluidity))



  fluidity_list <- list(fluidity=mean(all_comb$Fluidity),Standard_Deviation=stats::sd(all_comb$Fluidity),data=all_fluid,genome_fluidity=res)
  return(fluidity_list)}
