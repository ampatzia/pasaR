#'  grid plot
#'
#' This function outputs cluster tidyr::spread for genomes, genome participation per cluster,
#' gene participation per cluster and a brief summary allowing a quick exploration
#' of the available data
#'
#' @param Panmatrix Panmatrix produced by make_panmatrix functions
#' @param use_log Logical:Scale axis containing clusters with log
#' @export
#' @examples \dontrun{grid_plot(panm)}


grid_plot <- function(Panmatrix,use_log){

  if(missing(use_log)){use_log=TRUE}

  if(use_log==TRUE){
    a1 <- pm_plot(Panmatrix,use_log)+ggtitle(" Cluster spead for Genomes")+ylab("Clusters (log)")
    a2 <- cp_plot(Panmatrix,use_log)+ggtitle("Genome participation per Cluster")+xlab("Clusters (log)")
    a3 <- gp_plot(Panmatrix,use_log)+ggtitle("Gene participation per Cluster")+ylab("Clusters (log)")
    a4 <- mg_plot(Panmatrix,use_log)+ggtitle("Gene participation frequency")+ylab("Genes (log)")}else{

      a1 <- pm_plot(Panmatrix,use_log=FALSE)+ggtitle(" Cluster spead for Genomes")+ylab("Clusters")
      a2 <- cp_plot(Panmatrix,use_log=FALSE)+ggtitle("Genome participation per Cluster")+xlab("Clusters")
      a3 <- gp_plot(Panmatrix,use_log=FALSE)+ggtitle("Gene participation per Cluster")+ylab("Clusters")
      a4 <- mg_plot(Panmatrix,use_log)+ggtitle("Gene participation frequency")+ylab("Genes")
    }


  gridExtra::grid.arrange(a1, a2,a3,a4, ncol=2, top = "Panmatrix exploration Plots", padding = grid::unit(0.7, "line"))
}
