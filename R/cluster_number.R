#'  Cluster Analysis with fluidity
#'
#' This function outputs the proposed gene clustering based on fluidity distance
#'
#'
#' @param fluidity_list List produced by fluidity_all()
#' @param method  Method of hierarchical clustering, see ?stats::hclust()
#' @export
#' @examples \dontrun{cluster_number(fluidity_list)}
#'


cluster_number <- function(fluidity_list,method="ward.D"){

  Genome_1<-NULL
  Fluidity<-NULL
  dist_matrix <- tidyr::spread(fluidity_list$data,Genome_1,value=Fluidity,fill=0)
  dist_matrix <- dist_matrix[,-1]
  best_cluster <- data.frame(Clusters=NA,Index=NA,Value=NA)
  s1=NULL
  for(i in 2:10){
    k=stats:: cutree(stats::hclust(stats::as.dist(dist_matrix),method=method),i)
    s1[i]=fpc::cluster.stats(stats::as.dist(dist_matrix),k)$avg.silwidth
  }
  s1 <- s1[-1]
  best_cluster[1,] <- c(which.max(s1)+1,"Average Silhuette Width",round(max(s1, na.rm =TRUE),5))


  s2=NULL
  for(i in 2:10){
    k=stats:: cutree(stats::hclust(stats::as.dist(dist_matrix),method=method),i)
    s2[i]=fpc::cluster.stats(stats::as.dist(dist_matrix),k)$widestgap
  }
  s2<-s2[-1]
  best_cluster[2,]<-c(which.max(s2)+1,"Gap Statistic",round(max(s2, na.rm =TRUE),5))

  s3=NULL
  for(i in 2:10){
    k=stats:: cutree(stats::hclust(stats::as.dist(dist_matrix),method=method),i)
    s3[i]=fpc::cluster.stats(stats::as.dist(dist_matrix),k)$dunn
  }
  s3<-s3[-1]
  best_cluster[3,]<-c(which.max(s3)+1,"Dunn",round(max(s3, na.rm =TRUE),5))

  s4=NULL
  for(i in 2:10){
    k=stats:: cutree(stats::hclust(stats::as.dist(dist_matrix),method=method),i)
    s4[i]=fpc::cluster.stats(stats::as.dist(dist_matrix),k)$dunn
  }
  s4<-s4[-1]
  best_cluster[4,]<-c(which.min(s4)+1,"Entropy",round(min(s4, na.rm =TRUE),5))
  return(best_cluster)
}
