#'  Gene membership
#'
#' This function produces a genes - membership plot  of the pangenome, i.e. the gene-cluster participation
#' frequency in the pangenome.
#'
#' @param Panmatrix Panmatrix produced by make_panmatrix functions
#' @param collapsed Defaults to True, sum the number of genes for all groups with more than
#'                   three times the number of the genomes explored
#' @param use_log Logical,scales axis containing clusters with log, defaults to FALSE
#' @export
#' @examples \dontrun{mg_plot(Panmatrix)}
#' @importFrom dplyr n

mg_plot <- function(Panmatrix, collapsed ,use_log=TRUE){

  if (missing(collapsed)){
    collapsed = TRUE}

  Genes<-NULL
  n_memb <- colSums(Panmatrix)
  Cluster <- rep("Cluster",length(n_memb))

  sums <- data.frame(n_memb,Cluster)%>%
    dplyr::count(Cluster,n_memb)%>%dplyr::rename(.,Genes=n)


  if (collapsed == TRUE) {
    p_limit <- 3 * nrow(Panmatrix)
    y1 <- dplyr::filter(sums, n_memb >= p_limit)
    sums <- dplyr::filter(sums, n_memb < p_limit)
    y1 <- data.frame(n_memb= p_limit, n_memb = sum(y1$Genes))
    y3 <- dplyr::bind_rows(sums, y1)
  }

  if (use_log == TRUE) {sums$Genes <- log(sums$Genes)}


  p <- ggplot2::ggplot(sums, aes(x = n_memb, y = Genes))+ geom_point()+geom_line()+xlab("Number of Members")
  if (use_log == TRUE) {  p <- p + ylab("Genes (log)")}

  p
}
