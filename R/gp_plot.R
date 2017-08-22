#'  Gene participation plots
#'
#' This function produces a basic plot of the Gene paticipation per Cluster.
#'    This can be either a point plot or a bar plot.
#' @param Panmatrix Panmatrix produced by make_panmatrix functions
#' @param show_cluster Integer value  optional parameter allows to user to "zoom",
#'     ignoring clusters that have organism participation below it
#' @param plot_type Can be either line or bar
#' @param collapsed Logical, TRUE merges all clusters with more than 3 times the number of genomes examined
#' @param use_log Logical:Scale axis containing clusters with log
#' @export
#' @examples \dontrun{gp_plot(Panmatrix)}
#'


gp_plot <- function (Panmatrix, show_cluster, plot_type, collapsed=FALSE, use_log) {
  if (missing(show_cluster)) {show_cluster = 0}
  if (missing(plot_type)) {plot_type = "point"}
  if(missing(use_log)){use_log=TRUE}

  Genes<-NULL #prevent namespaces problems
  Cluster<-NULL

  levs <- 1:nrow(Panmatrix)
  y <- data.frame(Genes = colSums(Panmatrix), Cluster = seq(from = 1,
                                                            to = ncol(Panmatrix), by = 1))
  y <- as.data.frame(table(y$Genes))
  colnames(y)<-c("Genes", "Cluster")

  y$Genes <- as.numeric(as.character(y$Genes))
  if(collapsed==TRUE){
    p_limit <- 3*nrow(Panmatrix) #merge categories with more genes than 3* <number of genomes>
    y1 <- dplyr::filter(y, Genes>=p_limit)
    y <- dplyr::filter(y, Genes<p_limit)
    y1 <- data.frame(Genes=p_limit, Cluster=sum(y1$Cluster))
    y <- dplyr::bind_rows(y, y1)}

  if(use_log==TRUE){y$Cluster <- log(y$Cluster)}

  p <- ggplot2::ggplot(y, aes(x = Genes, y = Cluster))
  if(use_log==TRUE){p <- p+ylab("Cluster (log)")}

  if (plot_type == "bar") {
    p + geom_bar(stat = "identity")
  }
  else if (plot_type == "point") {
    p + geom_point()
  }
  else {
    message("This type of plot is not supported, use plot type point or bar.")
  }
}
