#'  Make panmatrix (MCL data)
#'
#' This function allows importing MCL output
#' @param file Disk Path to file
#' @export
#' @examples make_panmatrix()
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

#' @export
#' @examples make_panmatrix_fami()
#'


pm_plot<-function (object, show_cluster,plot_type)
{
  #show_cluster: optional parameter allows to user to "zoom", ignoring clusters that have organism participation
  #below it

  if(missing(show_cluster)){show_cluster=0}

  if(missing(plot_type)){plot_type="line"}


  object<-sapply(object,function(x) as.logical(x))%>%.[,!colSums(.)<show_cluster]

  levs <- 1:nrow(object)
  y <- as.data.frame(table(factor(colSums(object), levels = levs)))
  colnames(y)<-c("Genomes","Clusters")
  y$Genomes<-as.numeric(as.character(y$Genomes))
  y<-y[!y$Genomes<show_cluster,]
  p<-ggplot(y,aes(x=Genomes,y=Clusters))
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

#' @export
#' @examples cp_plot(panm)
#'


cp_plot<-function (object, show_cluster,plot_type)
{
  #show_cluster: optional parameter allows to user to "zoom", ignoring clusters that have organism participation
  #below it

  if(missing(show_cluster)){show_cluster=0}

  if(missing(plot_type)){plot_type="line"}


  object<-sapply(object,function(x) as.logical(x))%>%.[,!colSums(.)<show_cluster]

  levs <- 1:nrow(object)
  y <- data.frame(Genomes=colSums(object),Cluster=seq(from=1,to=ncol(object),by=1))

  colnames(y)<-c("Genomes","Clusters")
  y$Genomes<-as.numeric(as.character(y$Genomes))
  y<-y[!y$Genomes<show_cluster,]
  p<-ggplot(y,aes(y=Genomes,x=Clusters))+scale_y_continuous(breaks=seq(0,max(y$Genomes)+5,2))

  if(plot_type=="bar"){p+ geom_bar(stat="identity")}else if(plot_type=="point"){


    p+ geom_point()}else{
      message("This type of plot is not supported, use plot type point or bar.")
    }

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
#' @export
#' @examples gp_plot(panm)
#'


gp_plot<-function (object, show_cluster, plot_type,collapsed=FALSE) {
  if (missing(show_cluster)) {
    show_cluster = 0
  }
  if (missing(plot_type)) {
    plot_type = "point"
  }
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


  p <- ggplot(y, aes(x = Genes, y = Cluster))
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
#'    based on Heaps Law. This returns two estimated parameters, intercept and decay parameter a.
#'    If a<1 then the pangenome is considered to be open. This is an optimized version
#'    of the heaps() function from package micropan.
#'
#' @param Panmatrix Panmatrix produced by make_panmatrix functions
#' @param n.perm Number of permutations

#' @export
#' @examples pm_heaps(panm,100)
#'

pm_heaps<-function (panmatrix, n.perm)

{ pan.matrix<-sapply(panmatrix,function(x) as.logical(x))
ng<-nrow(panmatrix)
nmat <- matrix(0, nrow = (ng - 1), ncol = n.perm)

nmat<-replicate(n.perm,{

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
#' This function determines the closedness of the pangenome using a binomial mixture fit of the data
#'    based on Heaps Law. This returns pangenome and core size estimation. This is an optimized version
#'    of the binomixEstimate() function from package micropan.
#'
#' @param Panmatrix Panmatrix produced by make_panmatrix functions
#' @param K.range Range of model components to be tested

#' @export
#' @examples pm_binom(panm, K.range =2:8)
#'

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
#'    size along with variance and a 95% CI . This is a conservative estimator, ie. it is
#'    likely to produce a small estimation of the pangenome compared to other estimators.
#'    This function is an optimized and expanded version of chao() from package micropan.
#'
#'
#' @param Panmatrix Panmatrix produced by make_panmatrix functions
#' @export
#' @examples pm_chao(panm)
#'

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



#'  Genome Fluidity with sampling
#'
#' This function computes fluidity with sampling, as implemented on package micropan (optimized for speed).
#'   Fluidity takes values in [0,1], with 1 denoting no common genes.
#'   This metric was introduced in:
#'   "Genomic fluidity: an integrative view of gene diversity within microbial populations" (Kislyuk et al ,2011)
#'
#'
#'
#' @param Panmatrix Panmatrix produced by make_panmatrix functions
#' @param n.sim Number of simulations
#' @export
#' @examples pm_fluidity(panm, n.sim=100)
#'

pm_fluidity<-function (panm, n.sim = 10)
{
  panm<-sapply(as.data.frame(panm),function(x) as.logical(x))
  ng <- dim(panm)[1]
  flu <- rep(0, n.sim)
  for (i in 1:n.sim) {
    ii <- sample(ng, 2)
    g1 <- panm[ii[1], ]
    g2 <- panm[ii[2], ]
    flu[i] <- (sum(g1 > 0 & g2 == 0) + sum(g1 == 0 & g2 >
                                             0))/(sum(g1) + sum(g2))
  }
  flu.list <- list(Mean = mean(flu), Std = sd(flu))
  return(flu.list)
}


#'  Helper function: gtools::combination
#'
#' This helper function is produces all possible combinations for a set of numbers.
#'   This is an direct copy of the combination() function of package gtools, and is
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
#'   Fluidity takes values in [0,1], with 1 denoting no common genes.
#'   This metric was introduced in:
#'   "Genomic fluidity: an integrative view of gene diversity within microbial populations" (Kislyuk et al ,2011)
#'
#'
#'
#' @param Panmatrix Panmatrix produced by make_panmatrix functions
#' @param n.sim Number of simulations
#' @export
#' @examples pm_fluidity(panm)
#'

panm_fluidity_all<-function (panm){

  all_comb<-as.data.frame(gtools_comb(nrow(panm),2))
  panm <- sapply(panm, function(x) as.logical(x))

  fluid<-function(x){
    g1<-panm[x,] %>% .[,colSums(.)>0]
    flu<-(sum(g1[1,]==TRUE & g1[2,] == FALSE) + sum(g1[1,] == FALSE & g1[2,]==TRUE))/sum(colSums(g1))

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
#' This function outputs the genome names of a MVL clustering type input
#'
#'
#' @param Panmatrix Panmatrix produced by make_panmatrix functions
#' @param n.sim Number of simulations
#' @export
#' @examples organism_names_panmatrix(file)

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
  organism_names<-panm[,1]
  rm(panm)
  return(organism_names)

}


binomixMachine<-function (y, K, core.detect.prob = 1)
{
  n <- sum(y)
  G <- length(y)
  ctr <- list(maxit = 200*K, reltol = 1e-06)
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
