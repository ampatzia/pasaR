#'  Pangenome agglomerative hierarchical clustering based on fluidity
#'
#' This function performs agglomerative hierarchical clustering based on fluidity results
#' @param fluidity_list data produced from pm_fluidity_all
#' @param method method of clustering as used in stats::hclust()
#' @param genome_names optional file to label genomes as outputed by organism_names() and similar functions
#' @export
#' @examples \dontrun{pm_cluster(fluidity_result,"ward.D",genome_names)}

pm_cluster <- function(fluidity_list,method="ward.D",genome_names){

  Genome_1<-NULL
  Fluidity<-NULL

  dist_matrix <- tidyr::spread(fluidity_list$data, Genome_1, value = Fluidity,
                               fill = 0) #make matrix diagonal to convert to distance
  dist_matrix <- dist_matrix[, -1]
  clust_res <- stats::hclust(stats::as.dist(dist_matrix), method = method)
  if(class(genome_names)== "data.frame"){clust_res$labels <- genome_names$Organism}
  if(class(genome_names)== "character"){clust_res$labels <- genome_names}
  if(!missing(genome_names)){clust_res$labels <- paste0("gen_",seq(1,nrow(fluidity_list$genome_fluidity),1))}

  return(clust_res)
}
