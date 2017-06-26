#'  Make panmatrix (MCL data)
#'
#' This function allows importing MCL output
#' @param file Disk Path to file
#' @export
#' @examples make_panmatrix()
#' @note  MCL output as described in F. E. Psomopoulos, O. T. Vrousgou, and P. A. Mitkas, "Large-scale modular comparative genomics: the Grid approach [v1; not peer reviewed]," F1000research 2015, vol. 4(ISCB Com, iss. 377, p. 1, 2015. doi:10.7490/f1000research.1110127.1
#'        A. M. Kintsakis, F. E. Psomopoulos, and P. A. Mitkas, "Data-aware optimization of bioinformatics workflows in hybrid clouds," Journal of big data, vol. 3, iss. 20, pp. 1-26, 2016. doi:10.3389/fpls.2016.00554

#'
#'

make_panmatrix<-function(x){

make_base_df <- function(file){

  work_list <-scan(file=x,what="character,n=195,",sep=" ", allowEscapes = TRUE)%>%
    str_split_fixed(.," ", n = Inf) %>%
    sapply(.,stri_escape_unicode) %>% #Escapes all Unicode (not ASCII-printable) code points ie. single /
    sapply(., function(x) str_split_fixed(x,"[\\\\]+t|[^[:print:]]" , n = Inf)) %>%
    lapply(., function(x) str_split_fixed(x,"\\|", n=Inf))%>%
    lapply(., function(x) data.frame(x, stringsAsFactors=FALSE)) %>%
    lapply(., function(x){colnames(x)[1] <- "x1"; x}) %>%
    lapply(., function(x) separate(x,x1,into = c("Organism", "Protein", "Other"), sep="\\$"))%>%
    lapply(., function(x) x[,names(x) %in% c("Organism", "Protein") ])

  for (i in 1:length(work_list)){

    work_list[[i]]<-transform(work_list[[i]],Cluster=i)}


  #Make list to dataframe
  result_df<-bind_rows(work_list)
  rm(work_list)

  return(result_df)

}

cluster_composition<-function(x){
  result_df_cl1<-x %>% group_by(.,Cluster,Organism) %>%
    summarise(.,Proteins=length(Protein))
  return(result_df_cl1)}


  panm<-x %>% make_base_df(.) %>% cluster_composition(.) %>%spread(.,Cluster,Proteins,fill=0)
  #organism_names<-panm[,1]
  panm<-panm[,-1]
  return(panm)

}

#' Make panmatrix (fami MCL data)
#'
#' This function allows importing MCL output
#' @param path file Path to file

#' @export
#' @examples make_panmatrix_fami()
#'

make_panmatrix_fami<-function(file){
work_list <- read_delim(file,
                          "\t", escape_double = FALSE, col_names = FALSE,
                          trim_ws = TRUE)

split<-as.data.frame(str_split_fixed(work_list$X1," ",n=Inf))

work_list<-bind_cols(split,work_list[,-1])
work_list$V2<-as.character(work_list$V2)
work_list<-gather(work_list,V1)
work_list<-work_list[,c(1,3)]
colnames(work_list)<-c("Cluster","V1")
work_list<-work_list[complete.cases(work_list),]
work_list<-separate(work_list,V1,into=c("n1","n2","version","PID"),sep="-")
work_list$n1<-paste0(work_list$n1,"_",work_list$n2)
work_list$Cluster<-as.numeric(unlist(str_extract_all(work_list$Cluster,"\\(?[0-9,.]+\\)?")))

work_list<-work_list[,c(1,2,5)]
colnames(work_list)<-c("Cluster","Organism","Protein")

cluster_composition<-function(x){
  result_df_cl1<-x %>% group_by(.,Cluster,Organism) %>%
    summarise(.,Proteins=length(Protein))
  return(result_df_cl1)}



panm<-work_list %>% cluster_composition(.) %>%spread(.,Cluster,Proteins,fill=0)

org_names<-panm[,1]
panm<-panm[,-1]}

#' Make panmatrix (fami 2 MCL data)
#'
#' This function allows importing BLAST MCL output with default parameters
#' @param path file Path to file
#' @export
#' @examples make_panmatrix_fami2("path")
#'

make_panmatrix_fami2<-function (file){
  work_list <- read_delim(file, "\t", escape_double = FALSE,
                          col_names = FALSE, trim_ws = TRUE)

  work_list$V1<-(paste0("cluster",1:nrow(work_list)))
  colvals<-paste0("X",1:(ncol(work_list)-1))
  work_list <- gather(work_list,cluster,value,1:(ncol(work_list)-1))
  work_list <- work_list[, c(1, 3)]
  colnames(work_list) <- c("Cluster", "V1")
  work_list <- work_list[complete.cases(work_list), ]
  work_list <- separate(work_list, V1, into = c("n1", "n2",
                                                "version", "PID"), sep = "-")
  work_list$n1 <- paste0(work_list$n1, "_", work_list$n2)
  work_list$Cluster <- as.numeric(unlist(str_extract_all(work_list$Cluster,
                                                         "\\(?[0-9,.]+\\)?")))
  work_list <- work_list[, c(1, 2, 5)]
  colnames(work_list) <- c("Cluster", "Organism", "Protein")
  cluster_composition <- function(x) {
    result_df_cl1 <- x %>% group_by(., Cluster, Organism) %>%
      summarise(., Proteins = length(Protein))
    return(result_df_cl1)
  }
  panm <- work_list %>% cluster_composition(.) %>% spread(.,
                                                          Cluster, Proteins, fill = 0)
  org_names <- panm[, 1]
  panm <- panm[, -1]
}


#'  Panmatrix Summary
#'
#' This function produces a frequency table of Genome participation in clusters
#' @param Panmatrix  Panmatrix produced by make_panmatrix functions

#' @export
#' @examples make_panmatrix_fami()
#'

panm_summary<-function (object, ...)
{
  object<-sapply(object,function(x) as.logical(x))
  levs <- 1:nrow(object)
  y <- as.data.frame(table(factor(colSums(object), levels = levs)))
  colnames(y)<-c("Genomes","Clusters")
  return(y)
}


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
#' @examples pm_plot(panm)
#'


pm_plot<-function (object, show_cluster,plot_type,use_log)
{
  #show_cluster: optional parameter allows to user to "zoom", ignoring clusters that have organism participation
  #below it

  if(missing(show_cluster)){show_cluster=0}

  if(missing(plot_type)){plot_type="line"}
  if(missing(use_log)){use_log=TRUE}


  object<-sapply(object,function(x) as.logical(x))%>%.[,!colSums(.)<show_cluster]

  levs <- 1:nrow(object)
  y <- as.data.frame(table(factor(colSums(object), levels = levs)))
  colnames(y)<-c("Genomes","Clusters")
  y$Genomes<-as.numeric(as.character(y$Genomes))

  if(use_log==TRUE){y$Clusters<-log(y$Clusters)}

  y<-y[!y$Genomes<show_cluster,]
  p<-ggplot(y,aes(x=Genomes,y=Clusters))
  if(use_log==TRUE){p<-p+ylab("Cluster (log)")}
  if(plot_type=="bar"){p+ geom_bar(stat="identity")}else if(plot_type=="line"){


    p+ geom_line()+geom_point()}else{
      message("This type of plot is not supported, use plot type line or bar.")
    }

}




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
#' @examples cp_plot(panm)
#'


cp_plot<-function (object, show_cluster,plot_type,use_log)
{
  #s
  if(missing(show_cluster)){show_cluster=0}

  if(missing(plot_type)){plot_type="point"}
  if(missing(use_log)){use_log=TRUE}


  object<-sapply(object,function(x) as.logical(x))%>%.[,!colSums(.)<show_cluster]

  levs <- 1:nrow(object)
  y <- data.frame(Genomes=colSums(object),Cluster=seq(from=1,to=ncol(object),by=1))

  colnames(y)<-c("Genomes","Clusters")
  y$Genomes<-as.numeric(as.character(y$Genomes))

  if(use_log==TRUE){y$Clusters<-log(y$Clusters)}
  y<-y[!y$Genomes<show_cluster,]
  p<-ggplot(y,aes(y=Genomes,x=Clusters))+scale_y_continuous(breaks=seq(0,max(y$Genomes)+5,2))
  if(use_log==TRUE){p<-p+xlab("Cluster (log)")}

  if(plot_type=="bar"){p+ geom_bar(stat="identity")}else if(plot_type=="point"){


    p+ geom_point()}else{
      message("This type of plot is not supported, use plot type point or bar.")
    }
  #how_cluster: optional parameter allows to user to "zoom", ignoring clusters that have organism participation
  #below it


}



#'  Gene participation plots
#'
#' This function produces a basic plot of the Gene paticipation per Cluster.
#'    This can be either a point plot or a bar plot.
#' @param Panmatrix Panmatrix produced by make_panmatrix functions
#' @param show_cluster Integer value  optional parameter allows to user to "zoom",
#'     ignoring clusters that have organism participation below it
#' @param plot_type Can be either line or bar
#' @param plot_type Logical, TRUE merges all clusters with more than 3 times the number of genomes examined
#' @param use_log Logical:Scale axis containing clusters with log
#' @export
#' @examples gp_plot(panm)
#'


gp_plot<-function (object, show_cluster, plot_type,collapsed=FALSE,use_log) {
  if (missing(show_cluster)) {show_cluster = 0}
  if (missing(plot_type)) {plot_type = "point"}
  if(missing(use_log)){use_log=TRUE}

  levs <- 1:nrow(object)
  y <- data.frame(Genes = colSums(object), Cluster = seq(from = 1,
                                                         to = ncol(object), by = 1))
  y<-as.data.frame(table(y$Genes))
  colnames(y)<-c("Genes","Cluster")

  y$Genes<-as.numeric(as.character(y$Genes))
  if(collapsed==TRUE){
    p_limit<-3*nrow(object) #merge categories with more genes than 3* <number of genomes>
    y1<-filter(y,Genes>=p_limit)
    y<-filter(y,Genes<p_limit)
    y1<-data.frame(Genes=p_limit,Cluster=sum(y1$Cluster))
    y<-bind_rows(y,y1)}

  if(use_log==TRUE){y$Cluster<-log(y$Cluster)}

  p <- ggplot(y, aes(x = Genes, y = Cluster))
  if(use_log==TRUE){p<-p+ylab("Cluster (log)")}



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



#'  Heap's Law fitting
#'
#' This function determines the closedness of the pangenome using a curve fit of the data
#'    based on Heaps Law.
#'
#' @param Panmatrix Panmatrix produced by make_panmatrix functions
#' @param n_perm Number of permutations
#'
#' @details The regression fit returns two estimated parameters, intercept and decay parameter a.  If a<1 then the pangenome is considered to be open. This is an optimized version
#'    of the heaps() function from package micropan. The theoretical aspects are discussed in the canonical work of Tettelin et  al. (2008)
#'
#'
#'
#' @export
#' @examples pm_heaps(panm,100)
#'
#' @references  Tettelin, H., Riley, D., Cattuto, C., Medini, D. (2008). Comparative genomics: the bacterial pan-genome. Current Opinions in Microbiology, 12:472-477.

pm_heaps<-function (panmatrix, n_perm){
  if (missing(n_perm)) {n_perm = 100}
  pan.matrix<-sapply(panmatrix,function(x) as.logical(x))
ng<-nrow(panmatrix)
nmat <- matrix(0, nrow = (ng - 1), ncol = n_perm)

nmat<-replicate(n_perm,{

  cm <- apply(pan.matrix[sample(ng), ], 2, cumsum)
  rowSums((cm == 1)[2:ng, ] & (cm == 0)[1:(ng -  1), ])
})

nmat<-t(nmat)
colnames(nmat)<-c(2:(ncol(nmat)+1))
nmat<-gather(as.data.frame(nmat),genomes,genes)%>%transform(.,genomes=as.numeric(genomes))


p0 <- c(mean(nmat$genes[nmat$genomes == 2]), 1)


objectFun<-function (p, x, y)
{
  y.hat <- p[1] * x^(-p[2])
  J <- sqrt(sum((y - y.hat)^2))/length(x)
  return(J)
}

fit <- optim(p0, objectFun, gr = NULL, nmat$genomes, nmat$genes, method = "L-BFGS-B",
             lower = c(0, 0), upper = c(10000, 2))
p.hat <- fit$par
names(p.hat) <- c("Intercept", "alpha")
return(p.hat)
}


#'  Binomial mixture fitting
#'
#' This function determines the closedness of the pangenome using a binomial mixture fit of the data.
#'
#'
#' @param Panmatrix Panmatrix produced by make_panmatrix functions
#' @param K.range Range of model components to be tested
#'
#' @details The function returns pangenome and core size estimation along with the propabilities
#'  of the components.This is an optimized version of the binomixEstimate() function from package
#'  micropan  (Snipen & Hiland, 2015).
#'
#' @export
#' @examples pm_binom(panm, K.range =2:8)
#'
#' @references L. Snipen and K. H. Liland, "micropan: an R-package for microbial pan-genomics.," BMC bioinformatics, vol. 16, p. 79, 2015

pm_binom<-function (pan.matrix, K.range = 3:5, core.detect.prob = 1, verbose = TRUE)
{
  pan.matrix <- sapply(pan.matrix, function(x) as.logical(x))
  y <- table(factor(colSums(pan.matrix), levels = 1:dim(pan.matrix)[1]))
  bic.tab <- matrix(NA, nrow = length(K.range), ncol = 3)
  colnames(bic.tab) <- c("Core.size", "Pan.size", "BIC")
  rownames(bic.tab) <- paste(K.range, "components")
  mix.list <- vector("list", length(K.range))
  for (i in 1:length(K.range)) {
    lst <- binomixMachine(y, K.range[i], core.detect.prob)
    bic.tab[i, ] <- lst[[1]]
    mix.list[[i]] <- lst[[2]]
  }
  if (bic.tab[length(K.range), 3] == min(bic.tab[, 3]))
    warning("Minimum BIC at maximum K, increase upper limit of K.range")
  binomix <- list(BIC.table = bic.tab, Mix.list = mix.list)

  return(binomix)
}



pm_chao<-function (panm) {
  panm<-sapply(as.data.frame(panm),function(x) as.logical(x))
  y <- table(factor(colSums(panm), levels = 1:dim(panm)[1]))

  if (y[2] != 0){
    pan.size.biased <- round(sum(y) + y[1]^2/(2 * y[2]))
  }


  pan.size <- round(sum(y) + (y[1]*(y[1]-1))/(2 * (y[2]+1))) #is bias corrected
  fratio<-y[1]/y[2]
  chao.variance<-y[2]*(0.5*(fratio)^2 + fratio^3 + 0.25*fratio^4 ) #chao 1987

  C_var<- exp(1.96*sqrt(log(1+chao.variance/(pan.size-sum(y)))))

  l_bound<- sum(y) +(pan.size-sum(y))/C_var
  u_bound<- sum(y) +(pan.size-sum(y))*C_var

  if(is.null(pan.size.biased)) {pan.size.biased<-NA}

  results<-c(pan.size.biased,pan.size,chao.variance,l_bound,u_bound)

  names(results)<- c("Estimated pangenome size - Biased","Estimated pangenome size", " Estimator Variance","CI (95%) -Lower Bound","CI (95%)- Upper Bound")
  return(results)

}


#'  Chao lower bound estimator
#'
#' This function computes two versions of the lower bound Chao Estimator for the pangenome
#'    size.
#'
#'
#' @param Panmatrix Panmatrix produced by make_panmatrix functions
#'
#'@details The function prodused a biased & unbiased version of the metric along with variance
#'    and a 95% CI . This is a conservative estimator, ie. it is likely to produce a small
#'    estimation of the pangenome compared to other estimators.
#'    This function is an optimized and expanded version of chao() from package micropan.
#'
#'
#'
#'
#' @export
#' @examples pm_chao(panm)
#'
#' @references
#' Chao, A. (1987). Estimating the population size for capture-recapture data with unequal catchability. Biometrics 43, 783-791.
#' Gotelli, N.J. and Colwell, R.K., 2011. Estimating species richness. Biological diversity: frontiers in measurement and assessment, 12, pp.39-54.
#'
#'

pm_chao<-function (panm,biased=FALSE) {
  panm<-sapply(as.data.frame(panm),function(x) as.logical(x))
  y <- table(factor(colSums(panm), levels = 1:dim(panm)[1]))

  if (y[2] != 0){
    pan.size.biased <- round(sum(y) + y[1]^2/(2 * y[2]))
  }


  pan.size <- round(sum(y) + (y[1]*(y[1]-1))/(2 * (y[2]+1))) #is bias corrected
  fratio<-y[1]/y[2]
  chao.variance<-y[2]*(0.5*(fratio)^2 + fratio^3 + 0.25*fratio^4 ) #chao 1987

  C_var<- exp(1.96*sqrt(log(1+chao.variance/(pan.size-sum(y)))))

  l_bound<- sum(y) +(pan.size-sum(y))/C_var
  u_bound<- sum(y) +(pan.size-sum(y))*C_var

  if(is.null(pan.size.biased)) {pan.size.biased<-NA}
   if( biased==TRUE){
  results<-c(pan.size.biased,pan.size,chao.variance,l_bound,u_bound)

  names(results)<- c("Estimated pangenome size - Biased","Estimated pangenome size", " Estimator Variance","CI (95%) -Lower Bound","CI (95%)- Upper Bound")
  }else{
    results<-c(pan.size,chao.variance,l_bound,u_bound)

    names(results)<- c("Estimated pangenome size", " Estimator Variance","CI (95%) -Lower Bound","CI (95%)- Upper Bound")
  }



  return(results)

}



#'  Genome Fluidity with sampling
#'
#' This function computes fluidity with sampling, as implemented on package micropan (optimized for speed).
#'   Fluidity takes values in [0,1], with 1 denoting no common genes.

#'
#'
#'
#' @param Panmatrix Panmatrix produced by make_panmatrix functions
#' @param n.sim Number of simulations
#' @details Fluidity takes values in [0,1], with 1 denoting no common genes.This metric was introduced is called the Sorensen distance (Deza & Deza, 2009) and was first introduced in the
#' context of a pangenome analysis in (Kislyuk et al ,2011).
#' @export
#' @examples pm_fluidity(panm, n.sim=100)
#'
#' @references
#' A. O. Kislyuk, B. Haegeman, N. H. Bergman, and J. S. Weitz, "Genomic fluidity???: an integrative view of gene diversity within microbial populations," BMC genomics, pp. 12-32, 2011.
#' M. M. Deza and E. Deza, Encyclopedia of Distances. Springer, 2009.

pm_fluidity<-function (panm, n.sim = 10)
{
  panm<-sapply(as.data.frame(panm),function(x) as.logical(x))
  ng <- dim(panm)[1]
  flu <- rep(0, n.sim)
  for (i in 1:n.sim) {
    ii <- sample(ng, 2)
    g1 <- panm[ii[1], ]
    g2 <- panm[ii[2], ]
    flu[i] <- sum(abs(g1-g2))/(sum(g1) + sum(g2))
  }
  flu.list <- list(Mean = mean(flu), Std = sd(flu))
  return(flu.list)
}


#'  Helper function: gtools::combination
#'
#' This helper function is produces all possible combinations for a set of numbers.
#'   This is an direct copy of the combination() function from package gtools, and is
#'   used in function pm_fludity_all().
#'




gtools_comb<-function (n, r, v = 1:n, set = TRUE, repeats.allowed = FALSE)
{
  if (mode(n) != "numeric" || length(n) != 1 || n < 1 || (n%%1) !=
      0)
    stop("bad value of n")
  if (mode(r) != "numeric" || length(r) != 1 || r < 1 || (r%%1) !=
      0)
    stop("bad value of r")
  if (!is.atomic(v) || length(v) < n)
    stop("v is either non-atomic or too short")
  if ((r > n) & repeats.allowed == FALSE)
    stop("r > n and repeats.allowed=FALSE")
  if (set) {
    v <- unique(sort(v))
    if (length(v) < n)
      stop("too few different elements")
  }
  v0 <- vector(mode(v), 0)
  if (repeats.allowed)
    sub <- function(n, r, v) {
      if (r == 0)
        v0
      else if (r == 1)
        matrix(v, n, 1)
      else if (n == 1)
        matrix(v, 1, r)
      else rbind(cbind(v[1], Recall(n, r - 1, v)), Recall(n -
                                                            1, r, v[-1]))
    }
  else sub <- function(n, r, v) {
    if (r == 0)
      v0
    else if (r == 1)
      matrix(v, n, 1)
    else if (r == n)
      matrix(v, 1, n)
    else rbind(cbind(v[1], Recall(n - 1, r - 1, v[-1])),
               Recall(n - 1, r, v[-1]))
  }
  sub(n, r, v[1:n])
}



#'  Genome Fluidity
#'
#' This function computes fluidity without sampling.
#'
#'
#'
#' @param Panmatrix Panmatrix produced by make_panmatrix functions
#' @param n.sim Number of simulations
#' @export
#' @details  Fluidity takes values in [0,1], with 1 denoting no common genes.This metric was introduced is called the Sorensen distance (Deza & Deza, 2009) and was first introduced in the
#'   context of a pangenome analysis in (Kislyuk et al ,2011).
#'
#' @examples pm_fluidity_all(panm)
#'
#'  @references
#' A. O. Kislyuk, B. Haegeman, N. H. Bergman, and J. S. Weitz, "Genomic fluidity???: an integrative view of gene diversity within microbial populations," BMC genomics, pp. 12-32, 2011.
#' M. M. Deza and E. Deza, Encyclopedia of Distances. Springer, 2009.

pm_fluidity_all<-function (panm){

  all_comb<-as.data.frame(gtools_comb(nrow(panm),2))
  panm <- sapply(panm, function(x) as.logical(x))

  fluid<-function(x){
    g1<-panm[x,] %>% .[,colSums(.)>0]
    flu<-sum(abs(g1[1,]-g1[2,]))/sum(colSums(g1))

    return(flu)}

    all_comb$fluidity<-apply(all_comb,1,fluid)
    colnames(all_comb)<-c("Genome_1","Genome_2","Fluidity")



    test1<-split(all_comb,all_comb$Genome_1)
    test2<-split(all_comb,all_comb$Genome_2)
    dummy<-data.frame(Genome_1=as.numeric(NA),Genome_2=(NA),Fluidity=(NA))
    test1[[max(all_comb$Genome_2)]]<-dummy
    names(test1)[max(all_comb$Genome_2)]<-as.character(max(all_comb$Genome_2))


    test2$'1'<-dummy
    test1<-test1[order(as.numeric(names(test1)))]
    test2<-test2[order(as.numeric(names(test2)))]


    test2<-lapply(test2, function(x){x<-x[,c(2,1,3)]})
    test2<-lapply(test2,setNames,colnames(test1$`1`))

    tester_merge<-Map(rbind,test2, test1)


    all_fluid<-do.call(rbind,tester_merge)%>%arrange(.,Genome_1,Genome_2)
    all_fluid<-all_fluid[complete.cases(all_fluid),]
    res<-all_fluid%>% group_by(.,Genome_1)%>%summarise(.,Fluidity=mean(Fluidity))%>%arrange(.,desc(Fluidity))



    fluidity_list<-list(fluidity=mean(all_comb$Fluidity),Standard_Deviation=sd(all_comb$Fluidity),data=all_fluid,genome_fluidity=res)
    return(fluidity_list)}




#'  Cluster Analysis with fluidity
#'
#' This function outputs the proposed gene clustering based on fluidity distance
#'
#'
#' @param fluidity_list List produced by fluidity_all()
#' @param method  Method of hierarchical clustering, see ?hclust()
#' @export
#' @examples cluster_number(fluidity_list)
#'



cluster_number<-function(fluidity_list,method="ward.D"){

  dist_matrix<-spread(fluidity_list$data,Genome_1,value=Fluidity,fill=0)
  dist_matrix<-dist_matrix[,-1]
  best_cluster<-data.frame(Clusters=NA,Index=NA,Value=NA)
  s1=NULL
  for(i in 2:10){
    k=cutree(hclust(as.dist(dist_matrix),method=method),i)
    s1[i]=cluster.stats(as.dist(dist_matrix),k)$avg.silwidth
  }
  s1<-s1[-1]
  best_cluster[1,]<-c(which.max(s1)+1,"Average Silhuette Width",round(max(s1, na.rm =TRUE),5))


  s2=NULL
  for(i in 2:10){
    k=cutree(hclust(as.dist(dist_matrix),method=method),i)
    s2[i]=cluster.stats(as.dist(dist_matrix),k)$widestgap
  }
  s2<-s2[-1]
  best_cluster[2,]<-c(which.max(s2)+1,"Gap Statistic",round(max(s2, na.rm =TRUE),5))

  s3=NULL
  for(i in 2:10){
    k=cutree(hclust(as.dist(dist_matrix),method=method),i)
    s3[i]=cluster.stats(as.dist(dist_matrix),k)$dunn
  }
  s3<-s3[-1]
  best_cluster[3,]<-c(which.max(s3)+1,"Dunn",round(max(s3, na.rm =TRUE),5))

  s4=NULL
  for(i in 2:10){
    k=cutree(hclust(as.dist(dist_matrix),method=method),i)
    s4[i]=cluster.stats(as.dist(dist_matrix),k)$dunn
  }
  s4<-s4[-1]
  best_cluster[4,]<-c(which.min(s4)+1,"Entropy",round(min(s4, na.rm =TRUE),5))
  return(best_cluster)
}




#'  Organism Names - fami
#'
#' This function outputs the genome names of a fami clustering type input
#'
#'
#' @param Panmatrix Panmatrix produced by make_panmatrix functions
#' @param n.sim Number of simulations
#' @export
#' @examples organism_names_panmatrix_fami(file)
#'

organism_names_panmatrix_fami<-function (file) {
  work_list <- read_delim("~/Dataset#1/chlamydiae", "\t", escape_double = FALSE,
                          col_names = FALSE, trim_ws = TRUE)
  split <- as.data.frame(str_split_fixed(work_list$X1, " ",
                                         n = Inf))
  work_list <- bind_cols(split, work_list[, -1])
  work_list$V2 <- as.character(work_list$V2)
  work_list <- gather(work_list, V1)
  work_list <- work_list[, c(1, 3)]
  colnames(work_list) <- c("Cluster", "V1")
  work_list <- work_list[complete.cases(work_list), ]
  work_list <- separate(work_list, V1, into = c("n1", "n2",
                                                "version", "PID"), sep = "-")
  work_list$n1 <- paste0(work_list$n1, "_", work_list$n2)
  work_list$Cluster <- as.numeric(unlist(str_extract_all(work_list$Cluster,
                                                         "\\(?[0-9,.]+\\)?")))
  work_list <- work_list[, c(1, 2, 5)]
  colnames(work_list) <- c("Cluster", "Organism", "Protein")
  cluster_composition <- function(x) {
    result_df_cl1 <- x %>% group_by(., Cluster, Organism) %>%
      summarise(., Proteins = length(Protein))
    return(result_df_cl1)
  }
  panm <- work_list %>% cluster_composition(.) %>% spread(.,
                                                          Cluster, Proteins, fill = 0)
  org_names <- panm[, 1]
  return(org_names)
}


#'  Organism Names - MCL
#'
#' This function outputs the genome names of a MCL clustering type input
#'
#'
#' @param Panmatrix Panmatrix produced by make_panmatrix functions
#' @export
#' @examples organism_names_panmatrix("path")

org_names<-function(x){

  make_base_df <- function(file){

    work_list <-scan(file=x,what="character,n=195,",sep=" ", allowEscapes = TRUE)%>%
      str_split_fixed(.," ", n = Inf) %>%
      sapply(.,stri_escape_unicode) %>% #Escapes all Unicode (not ASCII-printable) code points ie. single /
      sapply(., function(x) str_split_fixed(x,"[\\\\]+t|[^[:print:]]" , n = Inf)) %>%
      lapply(., function(x) str_split_fixed(x,"\\|", n=Inf))%>%
      lapply(., function(x) data.frame(x, stringsAsFactors=FALSE)) %>%
      lapply(., function(x){colnames(x)[1] <- "x1"; x}) %>%
      lapply(., function(x) separate(x,x1,into = c("Organism", "Protein", "Other"), sep="\\$"))%>%
      lapply(., function(x) x[,names(x) %in% c("Organism", "Protein") ])

    for (i in 1:length(work_list)){

      work_list[[i]]<-transform(work_list[[i]],Cluster=i)}


    #Make list to dataframe
    result_df<-bind_rows(work_list)
    rm(work_list)

    return(result_df)

  }

  cluster_composition<-function(x){
    result_df_cl1<-x %>% group_by(.,Cluster,Organism) %>%
      summarise(.,Proteins=length(Protein))
    return(result_df_cl1)}

  panm<-x %>% make_base_df(.) %>% cluster_composition(.) %>%spread(.,Cluster,Proteins,fill=0)
  organism_names<-panm[,1]
  rm(panm)
  return(organism_names)

}


#'  Organism Names - MCL (fami2)
#'
#' This function outputs the genome names of a MCL clustering type input
#'
#'
#' @param Panmatrix Panmatrix produced by make_panmatrix functions
#' @export
#' @examples organism_names_fami2("path")


org_names_fami2<-function (file){
  work_list <- read_delim(file, "\t", escape_double = FALSE,
                          col_names = FALSE, trim_ws = TRUE)

  work_list$V1<-(paste0("cluster",1:nrow(work_list)))
  colvals<-paste0("X",1:(ncol(work_list)-1))
  work_list <- gather(work_list,cluster,value,1:(ncol(work_list)-1))
  work_list <- work_list[, c(1, 3)]
  colnames(work_list) <- c("Cluster", "V1")
  work_list <- work_list[complete.cases(work_list), ]
  work_list <- separate(work_list, V1, into = c("n1", "n2",
                                                "version", "PID"), sep = "-")
  work_list$n1 <- paste0(work_list$n1, "_", work_list$n2)
  work_list$Cluster <- as.numeric(unlist(str_extract_all(work_list$Cluster,
                                                         "\\(?[0-9,.]+\\)?")))
  work_list <- work_list[, c(1, 2, 5)]
  colnames(work_list) <- c("Cluster", "Organism", "Protein")
  cluster_composition <- function(x) {
    result_df_cl1 <- x %>% group_by(., Cluster, Organism) %>%
      summarise(., Proteins = length(Protein))
    return(result_df_cl1)
  }
  panm <- work_list %>% cluster_composition(.) %>% spread(.,
                                                          Cluster, Proteins, fill = 0)
  org_names <- panm[, 1]
  return(org_names)
}

#'  Binomix machine
#'
#' This function is a helper borrowed from package micropan to be used to compute
#'  binomial mixture pangenome models
#'


binomixMachine<-function (y, K, core.detect.prob = 1)
{
  n <- sum(y)
  G <- length(y)
  ctr <- list(maxit = 200*K, reltol = 1e-8)
  np <- K - 1
  pmix0 <- rep(1, np)/K
  pdet0 <- (1:np)/(np + 1)
  p.initial <- c(pmix0, pdet0)
  A <- rbind(c(rep(1, np), rep(0, np)), c(rep(-1, np), rep(0,
                                                           np)), diag(np + np), -1 * diag(np + np))
  b <- c(0, -1, rep(0, np + np), rep(-1, np + np))
  est <- constrOptim(theta = p.initial, f = negTruncLogLike,
                     grad = NULL, method = "Nelder-Mead", control = ctr, ui = A,
                     ci = b, y = y, core.p = core.detect.prob)
  estimates <- numeric(3)
  names(estimates) <- c("Core.size", "Pan.size", "BIC")
  estimates[3] <- 2 * est$value + log(n) * (np + K)
  p.mix <- c(1 - sum(est$par[1:np]), est$par[1:np])
  p.det <- c(core.detect.prob, est$par[(np + 1):length(est$par)])
  ixx <- order(p.det)
  p.det <- p.det[ixx]
  p.mix <- p.mix[ixx]
  theta_0 <- choose(G, 0) * sum(p.mix * (1 - p.det)^G)
  y_0 <- n * theta_0/(1 - theta_0)
  estimates[2] <- n + round(y_0)
  ixx <- which(p.det >= core.detect.prob)
  estimates[1] <- round(estimates[2] * sum(p.mix[ixx]))
  mixmod <- matrix(c(p.det, p.mix), nrow = 2, byrow = T)
  rownames(mixmod) <- c("Detection.prob", "Mixing.prop")
  colnames(mixmod) <- paste("Comp_", 1:K, sep = "")
  return(list(estimates, mixmod))
}

#'   Truncated log likelihood
#'
#' This function is a helper borrowed from package micropan to be used to compute
#'  truncated log likelihood
#'

negTruncLogLike<-function (p, y, core.p)
{
  np <- length(p)/2
  p.det <- c(core.p, p[(np + 1):length(p)])
  p.mix <- c(1 - sum(p[1:np]), p[1:np])
  G <- length(y)
  K <- length(p.mix)
  n <- sum(y)
  theta_0 <- choose(G, 0) * sum(p.mix * (1 - p.det)^G)
  L <- -n * log(1 - theta_0)
  for (g in 1:G) {
    theta_g <- choose(G, g) * sum(p.mix * p.det^g * (1 -
                                                       p.det)^(G - g))
    L <- L + y[g] * log(theta_g)
  }
  return(-L)
}



#'  grid plot
#'
#' This function outputs cluster spread for genomes, genome participation per cluster,
#' gene participation per cluster and a brief summary allowing a quick exploration
#' of the available data
#'
#' @param Panmatrix Panmatrix produced by make_panmatrix functions
#' @param use_log Logical:Scale axis containing clusters with log
#' @export
#' @examples grid_plot(panm)


grid_plot<-function(panm,use_log){

  if(missing(use_log)){use_log=TRUE}

  if(use_log==TRUE){
  a1<-pm_plot(panm,use_log)+ggtitle(" Cluster spead for Genomes")+ylab("Clusters (log)")
  a2<-cp_plot(panm,use_log)+ggtitle("Genome participation per Cluster")+xlab("Clusters (log)")
  a3<-gp_plot(panm,use_log)+ggtitle("Gene participation per Cluster")+ylab("Clusters (log)")}else{

        a1<-pm_plot(panm,use_log=FALSE)+ggtitle(" Cluster spead for Genomes")+ylab("Clusters")
        a2<-cp_plot(panm,use_log=FALSE)+ggtitle("Genome participation per Cluster")+xlab("Clusters")
        a3<-gp_plot(panm,use_log=FALSE)+ggtitle("Gene participation per Cluster")+ylab("Clusters")
           }

  a4<-panm_summary(panm)
          if(nrow(a4)>14){
            cut_seq<-seq(0,nrow(panm)+1,8)
            cut_seq[length(cut_seq)]<-nrow(panm)
            a4$Genomes<-as.numeric(as.factor(a4$Genomes))
            a4$Genomes<-cut(a4$Genomes,cut_seq)
            a4<-a4%>%group_by(Genomes) %>%
            summarise(Clusters = sum(Clusters))}

  a4<-tableGrob(a4,theme=ttheme_minimal(base_size = 9))

  grid.arrange(a1, a2,a3,a4, ncol=2, top = "Panmatrix exploration Plots", padding = unit(0.7, "line"))
  }

#'  Pangenome agglomerative hierarchical clustering based on fluidity
#'
#' This function
#' @param fluidity_list data produced from pm_fluidity_all
#' @param method method of clustering as used in hclust()
#' @param genome_names optional file to label genomes as outputed by organism_names() and similar functions
#' @export
#' @examples pm_cluster(fluidity_result,"ward.D",genome_names)

pm_cluster<-function(fluidity_list,method="ward.D",genome_names){

  dist_matrix <- spread(fluidity_list$data, Genome_1, value = Fluidity,
                        fill = 0) #make matrix diagonal to convert to distance
  dist_matrix <- dist_matrix[, -1]
  clust_res<-hclust(as.dist(dist_matrix), method = method)
  if(!missing(genome_names)){clust_res$labels<-genome_names$Organism}
  return(clust_res)
}
