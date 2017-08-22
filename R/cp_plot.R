#'  Cluster participation plots
#'
#' This function produces a basic plot of the pangenome.
#'    This can be either a point plot or a bar plot.
#' @param Panmatrix Panmatrix produced by make_panmatrix functions
#' @param show_cluster Integer value  optional parameter allows to user to "zoom",
#'     ignoring clusters that have organism participation below it
#' @param plot_type Can be either line or bar
#' @param use_log Logical:Scale axis containing clusters with log
#' @export
#' @examples \dontrun{cp_plot(Panmatrix)}
#'


cp_plot <- function (Panmatrix, show_cluster,plot_type,use_log)
{
  if(missing(show_cluster)){show_cluster=0}

  if(missing(plot_type)){plot_type="point"}
  if(missing(use_log)){use_log=TRUE}
  Genomes<-NULL
  Clusters<-NULL

  object <- sapply(Panmatrix, function(x) as.logical(x))%>%.[, !colSums(.)<show_cluster]

  levs <- 1:nrow(object)
  y <- data.frame(Genomes=colSums(object), Cluster=seq(from=1, to=ncol(object), by=1))

  colnames(y) <- c("Genomes", "Clusters")
  y$Genomes <- as.numeric(as.character(y$Genomes))

  if(use_log==TRUE){y$Clusters <- log(y$Clusters)}
  y <- y[!y$Genomes<show_cluster,]
  p <- ggplot2::ggplot(y, aes(y=Genomes,x=Clusters))+scale_y_continuous(breaks=seq(0, max(y$Genomes)+5, 2))
  if(use_log==TRUE){p<-p+xlab("Cluster (log)")}

  if(plot_type=="bar"){p+ geom_bar(stat="identity")}else if(plot_type=="point"){


    p+ geom_point()}else{
      message("This type of plot is not supported, use plot type point or bar.")
    }
  #how_cluster: optional parameter allows to user to "zoom", ignoring clusters that have organism participation
  #below it


}
