#'  Panmatrix Plots
#'
#' This function produces a basic plot of the pangenome.
#'    This can be either a line plot or a bar plot.
#' @param Panmatrix Panmatrix produced by make_panmatrix functions
#' @param show_cluster Integer value  optional parameter allows to user to "zoom",
#'     ignoring clusters that have organism participation below it
#' @param plot_type Can be either line or bar
#' @param use_log Logical:Scale axis containing clusters with log
#' @export
#' @examples \dontrun{pm_plot(Panmatrix)}
#'


pm_plot <- function (Panmatrix, show_cluster,plot_type, use_log)
{
  #show_cluster: optional parameter allows to user to "zoom", ignoring clusters that have organism participation
  #below it

  if(missing(show_cluster)){show_cluster=0}

  if(missing(plot_type)){plot_type="line"}
  if(missing(use_log)){use_log=TRUE}

  Genomes<-NULL #prevent namespace problems
  Clusters<-NULL

  object <- sapply(Panmatrix, function(x) as.logical(x))%>%.[, !colSums(.)<show_cluster]

  levs <- 1:nrow(object)
  y <- as.data.frame(table(factor(colSums(object), levels = levs)))
  colnames(y) <- c("Genomes","Clusters")
  y$Genomes <- as.numeric(as.character(y$Genomes))

  if(use_log==TRUE){y$Clusters<-log(y$Clusters)}

  y <- y[!y$Genomes<show_cluster,]
  p <- ggplot2::ggplot(y,aes(x=Genomes, y=Clusters))
  if(use_log==TRUE){p <- p+ylab("Cluster (log)")}

  if(plot_type=="bar"){p+ geom_bar(stat="identity")}else if(plot_type=="line"){


    p+ geom_line()+geom_point()}else{
      message("This type of plot is not supported, use plot type line or bar.")
    }

}
