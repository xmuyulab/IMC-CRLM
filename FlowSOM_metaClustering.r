
#' PLEASE NOTE:
#' This file is one part of FlowSOM R package(4_metaClustering.R) , Please don't forget to cite the original FlowSOM article,if it is used.
#' I made some modifications as follows:
#' 1) add PhenoGraph as a new metacluster method
#' 2) add k related parameters,such as kmax(used to named by "max"),kstep,kmin
#' 3) use blue colored cross mark the elbow point in the curve
#' 4) add the function of metaClusting, which could did metaclustering,elbow_test and tsne_visualiaztion together 
#' Xinlei Chen, 2020-2-4


#' MetaClustering
#'
#' Cluster data with automatic number of cluster determination for 
#' several algorithms
#'
#' @param data   Matrix containing the data to cluster
#' @param method Clustering method to use
#' @param kmax    Maximum number of clusters to try out
#' @param elbow_test  Bool value, determine whether perform elbow_test
#' @param ...    Extra parameters to pass along
#' 
#' @return Numeric array indicating cluster for each datapoint
#' @seealso   \code{\link{metaClustering_consensus}}
#'
#' @examples
#'    # Read from file, build self-organizing map and minimal spanning tree
#'    fileName <- system.file("extdata", "68983.fcs", package="FlowSOM")
#'    flowSOM.res <- ReadInput(fileName, compensate=TRUE,transform=TRUE,
#'                             scale=TRUE)
#'    flowSOM.res <- BuildSOM(flowSOM.res,colsToUse=c(9,12,14:18))
#'    flowSOM.res <- BuildMST(flowSOM.res)
#'    
#'    # Apply metaclustering
#'    metacl <- MetaClustering(flowSOM.res$map$codes,
#'                             "metaClustering_consensus",
#'                             max=10)
#'    
#'    # Get metaclustering per cell
#'    flowSOM.clustering <- metacl[flowSOM.res$map$mapping[,1]]    
#'
#' @export
MetaClustering <- function(data,method,elbow_test=T,k_value=NULL,kmax=NULL,kstep=NULL,
                           kmin=NULL,plot = T,...){
  if(is.null(k_value)&(elbow_test==F)){
    elbow_test=T
    cat("k_value is not specified, start elbow test...\n")
  }
  res<-k_value  
   if(elbow_test==T){
         res <- DetermineNumberOfClusters(data,kmax,method,...)
   }
    method <- get(method)
    cat(paste0("k_value is set to ",res,"\n"))
    method(data,k=res)
}



DetermineNumberOfClusters <- function(data,
                                      kmax=NULL,
                                      method,
                                      plot=TRUE,
                                      smooth=0.2,
                                      kstep=NULL,
                                      kmin=NULL,
                                      ...){
  

  if(method=="metaClustering_PhenoGraph"){

    if(is.null(kmin))  kmin=5
    if(is.null(kstep)) kstep=5
    if(is.null(kmax))   kmax=90
    kseq<-seq(from=kmin,to=kmax,by=kstep)
    print('xxxxxxxx')
    print(kseq)
    res <- rep(0,length(kseq))
    nclus<-rep(0,length(kseq))
    
    for(i in 1:length(kseq)){
      #cat(paste0("PhenoGraph k is set to"),i,"...")
      #c <- as.numeric(membership(Rphenograph(data,k=kseq[i]))) xu
      c <- as.numeric(Rphenograph(data,k=kseq[i])[[2]]$membership)
      #cat(paste0(" Clustering Finished, Get", c, "clusters\n"))
      nclus[i]<-max(c)
      #res[i] <- SSE(data,c)
    }
    
    #Elbowpoint1 <-findElbow2(res,nclus)
    Elbowpoint2 <-findElbow2(nclus,kseq)
    
    #output plot with the elbowpoint labeled
    xadjust<-(max(kseq)-min(kseq))*0.07
    yadjust<-(max(nclus)-min(nclus))*0.07
    
    if(plot){# par(mfrow=c(1,2))
                   # plot(nclus, res, type="b", xlab="Number of Clusters",ylab="Within groups sum of squares")
                   # points(nclus[Elbowpoint1],res[Elbowpoint1],pch=3,col="blue",cex=2)
                   # text(nclus[Elbowpoint1]*0.8,res[Elbowpoint1],labels = paste0("x=",Elbowpoint1))
                          
                    plot(kseq, nclus, type="b", xlab="PhenGraph_k",ylab="Number of Clusters",title=paste0(method," Elbow Test"))
                    points(kseq[Elbowpoint2],nclus[Elbowpoint2],pch=3,col="blue",cex=2)
                    text(kseq[Elbowpoint2]+xadjust,nclus[Elbowpoint2]+yadjust,labels = paste0("x=",kseq[Elbowpoint2]))
 
      
    }
    cat(paste0("\nFind elbow point: PhenoGraph k =",kseq[Elbowpoint2],"\n" ))
    return(kseq[Elbowpoint2])
  }
  
  
  # if(method == "metaClustering_consensus"){
  #       if(is.null(kmin))  kmin=2
  #       if(is.null(kstep)) kstep=1
  #       if(is.null(max))   max=20
  #       kseq<-seq(from=kmin,to=max,by=kstep)
  #       res <- rep(0,length(kseq))
  #       results <- consensus(data,max,...)
  #       res <- rep(0,max)
  #       res[1] <- SSE(data,rep(1,nrow(data)))
  #       
  #       for(i in seq(from=kmin,to=max,by=kstep)){
  #           c <- results[[i]]$consensusClass
  #           res[i] <- SSE(data,c)
  #       }
  
  if(method ==    "metaClustering_consensus"){
        methodname<-method
    
        if(is.null(kmax))   kmax=20
        if(!is.null(kmin)|!is.null(kstep)) {cat("kmin and kstep can only set to 1 when clustering with metaClustering_consensus ")}
        kseq<-seq(from=1,to=kmax,by=1)
        results <- consensus(data,kmax,...)
        res <- rep(0,kmax)
        res[1] <- SSE(data,rep(1,nrow(data)))
        for(i in 2:kmax){
          c <- results[[i]]$consensusClass
          res[i] <- SSE(data,c)
        }
  }else
    {
        if(is.null(kmin))  kmin=1
        if(is.null(kstep)) kstep=1
        if(is.null(kmax))   kmax=20
        kseq<-seq(from=kmin,to=kmax,by=kstep)
        methodname<-method
        method <- get(method)
        res <- rep(0,length(kseq))
        for(i in 1:length(kseq)){
          c <- method(data, k=kseq[i])#修改一处bugc <- method(data, k=i,...)
          res[i] <- SSE(data,c)
        }
    }
  
    for(i in 2:(length(kseq)-1)){
        res[i] <- (1-smooth)*res[i]+(smooth/2)*res[i-1]+(smooth/2)*res[i+1]
    }
    
    Elbowpoint <-findElbow2(res,kseq)

    #output plot with the elbowpoint labeled
    xadjust<-(max(kseq)-min(kseq))*0.07
    yadjust<-(max(res)-min(res))*0.07
    
    
    if(plot){ plot(kseq, res, type="b", xlab="Number of Clusters", 
                    ylab="Within groups sum of squares",
                    title=paste0(methodname," Elbow Test"))
              points(kseq[Elbowpoint],res[Elbowpoint],pch=3,col="blue",cex=2)
              text(kseq[Elbowpoint]+xadjust,res[Elbowpoint]+yadjust,labels = paste0("x=",kseq[Elbowpoint]))
              
    }
    cat(paste0("Find elbow point: Cluster Number=",kseq[Elbowpoint],"\n" ))
    
    return(kseq[Elbowpoint])
}

findElbow <- function(data){
    n <- length(data)    
    data <- as.data.frame(cbind(1:n,data))
    colnames(data) <- c("X","Y")
    
    min_r <- Inf
    optimal <- 1
    for(i in 2:(n-1)){
        f1 <- stats::lm(Y~X,data[1:(i-1),])
        f2 <- stats::lm(Y~X,data[i:n,])
        r <- sum(abs(c(f1$residuals,f2$residuals)))
        if(r < min_r){
            min_r <- r
            optimal <-i
        }
    }
    optimal
}

findElbow2 <- function(data,kseq){
  n <- length(data)    
  data <- as.data.frame(cbind(kseq,data))
  colnames(data) <- c("X","Y")
  
  min_r <- Inf
  optimal <- 1
  for(i in 2:(n-1)){
    f1 <- stats::lm(Y~X,data[1:(i-1),])
    f2 <- stats::lm(Y~X,data[i:n,])
    r <- sum(abs(c(f1$residuals,f2$residuals)))
    if(r < min_r){
      min_r <- r
      optimal <-i
    }
  }
  optimal
}





#' MetaClustering
#' 
#' Cluster data using hierarchical consensus clustering with k clusters
#'
#' @param data Matrix containing the data to cluster
#' @param k    Number of clusters
#' @param seed Seed to pass to consensusClusterPlus
#' 
#' @return  Numeric array indicating cluster for each datapoint
#' @seealso \code{\link{MetaClustering}}
#' @examples
#'    # Read from file, build self-organizing map and minimal spanning tree
#'    fileName <- system.file("extdata", "68983.fcs", package="FlowSOM")
#'    flowSOM.res <- ReadInput(fileName, compensate=TRUE,transform=TRUE,
#'                             scale=TRUE)
#'    flowSOM.res <- BuildSOM(flowSOM.res,colsToUse=c(9,12,14:18))
#'    flowSOM.res <- BuildMST(flowSOM.res)
#'    
#'    # Apply consensus metaclustering
#'    metacl <- metaClustering_consensus(flowSOM.res$map$codes,k=10)    
#'
#' @export
metaClustering_consensus <- function(data, k=7,seed=NULL){
    results <- suppressMessages(ConsensusClusterPlus::ConsensusClusterPlus(
                                t(data),
                                maxK=k, reps=100, pItem=0.9, pFeature=1, 
                                title=tempdir(), plot="pdf", verbose=FALSE,
                                clusterAlg="hc", # "hc","km","kmdist","pam"
                                distance="euclidean" ,
    #"euclidean","pearson","spearman","binary","maximum","canberra","minkowski"
                                seed=seed
    ))
    
    results[[k]]$consensusClass
}

consensus <- function(data,max,...){
    results <- suppressMessages(ConsensusClusterPlus::ConsensusClusterPlus(
                                t(data),
                                maxK=max, reps=100, pItem=0.9, pFeature=1,
                                title=tempdir(), plot="pdf", verbose=FALSE,
                                clusterAlg="hc", # "hc","km","kmdist","pam"
                                distance="euclidean" 
    #"euclidean","pearson","spearman","binary","maximum","canberra","minkowski"
    ))
}

metaClustering_hclust <- function(data, k=7){
    d <- stats::dist(data, method = "minkowski")
    fit <- stats::hclust(d, method="ward.D2")
    stats::cutree(fit, k=k)
}

metaClustering_kmeans <- function(data, k=7){
    stats::kmeans(data, centers=k)$cluster
}

metaClustering_som <- function(data, k=7){
    s <- SOM(data,xdim=k,ydim=1,rlen=100)
    s$unit.classif
}


metaClustering_PhenoGraph<-function(data,k=30){
  #as.numeric(membership(Rphenograph(data,k=k))) xu
  as.numeric(Rphenograph(data,k=k)[[2]]$membership)
}



SSE <- function(data,clustering){
    if(class(clustering)!= "numeric")
        clustering <- as.numeric(as.factor(clustering))
    c_wss <- 0
    for(j in seq_along(clustering)){
        if(sum(clustering==j) > 1){
            c_wss <- c_wss + (nrow(data[clustering==j,,drop=FALSE])-1)*
                        sum(apply(data[clustering==j,,drop=FALSE],2,stats::var))
        }
    }
    c_wss
}



#' F measure
#' 
#' Compute the F measure between two clustering results
#'
#' @param realClusters Array containing real cluster labels for each sample
#' @param predictedClusters Array containing predicted cluster labels for each
#'                          sample
#' @param silent    Logical, if FALSE (default), print some information about 
#'                  precision and recall
#' 
#' @return  F measure score
#' @examples
#' # Generate some random data as an example
#' realClusters <- sample(1:5,100,replace = TRUE)
#' predictedClusters <- sample(1:6, 100, replace = TRUE)
#' 
#' # Calculate the FMeasure
#' FMeasure(realClusters,predictedClusters)
#' @export
FMeasure <- function(realClusters, predictedClusters,silent=FALSE){
    if (sum(predictedClusters)==0)
        return(0);
    a <- table(realClusters, predictedClusters);
    p <- t(apply(a,1,function(x)x/colSums(a)))
    r <- apply(a,2,function(x)x/rowSums(a))
    if(!silent) message("Precision: ",
                sum(apply(p,1,max) * (rowSums(a)/sum(a))),
                "\nRecall: ",sum(apply(r,1,max) * (rowSums(a)/sum(a))),"\n")
    f <- 2*r*p / (r+p)
    f[is.na(f)] <- 0
    sum(apply(f,1,max) * (rowSums(a)/sum(a)))
}

#' MetaclusterMFIs
#' 
#' Compute the median fluorescence intensities for the metaclusters
#'
#' @param fsom Result of calling the FlowSOM function
#' @return  Metacluster MFIs
#' @examples
#' fileName <- system.file("extdata", "68983.fcs", package="FlowSOM")
#' ff <- flowCore::read.FCS(fileName)
#' ff <- flowCore::compensate(ff,ff@@description$SPILL)
#' ff <- flowCore::transform(ff,
#'          flowCore::transformList(colnames(ff@@description$SPILL),
#'                                 flowCore::logicleTransform()))
#' flowSOM.res <- FlowSOM(ff,scale=TRUE,colsToUse=c(9,12,14:18),maxMeta=10)
#' mfis <- MetaclusterMFIs(flowSOM.res)
#' @export
MetaclusterMFIs <- function(fsom){
  MFIs <- t(sapply(seq_along(levels(fsom$metaclustering)), 
                  function(i) {
                    apply(subset(fsom$FlowSOM$data, 
                                 fsom$metaclustering[
                                   fsom$FlowSOM$map$mapping[,1]] == i),
                          2,
                          stats::median)
                  }))
  rownames(MFIs) <- seq_len(nrow(MFIs))
  return(MFIs)
}

#' MetaclusterCVs
#' 
#' Compute the coefficient of variation for the metaclusters
#'
#' @param fsom Result of calling the FlowSOM function
#' @return  Metacluster CVs
#' @examples
#' fileName <- system.file("extdata", "68983.fcs", package="FlowSOM")
#' ff <- flowCore::read.FCS(fileName)
#' ff <- flowCore::compensate(ff,ff@@description$SPILL)
#' ff <- flowCore::transform(ff,
#'          flowCore::transformList(colnames(ff@@description$SPILL),
#'                                 flowCore::logicleTransform()))
#' flowSOM.res <- FlowSOM(ff,scale=TRUE,colsToUse=c(9,12,14:18), nClus=10)
#' cvs <- MetaclusterCVs(flowSOM.res)
#' @export
MetaclusterCVs <- function(fsom){
  CVs <- t(sapply(seq_along(levels(fsom$metaclustering)), 
                  function(i) {
                    apply(subset(fsom$FlowSOM$data, 
                                 fsom$metaclustering[
                                   fsom$FlowSOM$map$mapping[,1]] == i),
                          2,
                          function(y){
                            if(length(y) > 0 && mean(y) != 0){
                              stats::sd(y)/mean(y)
                            } else {
                              NA
                            }})
                  }))
  return(CVs)
}



#' MetaClustering, elbow_test and tsne visualisation
#'
#' 
#' @param indataframe            event expression matrix or dataframe， containing one clustering channel
#' @param clustername            string,the column name of clustering channel
#' @param metaClustering_method  string,the method of metaclustering to use. could be set to "metaClustering_consensus",
#'                               "metaClustering_hclust","metaClustering_kmeans","metaClustering_PhenoGraph"
#' @param usecol                 numeric, the numeric columns id of indataframe, data in selected columns would be used in metaclustering. 
#' @param elbow_test             TRUE/FALSE, determine whether elbow_test will be carried out before formal metaclustering.
#' @param k_value                numeric, the k value which is used in the metaclustering_method, equal to the number of generated clusters except for PhenoGraph. 
#' @param view_tsne              TRUE/FALSE, determine whether tsne visualisation is performed, use TRUE as default.                               
#' @param seed,perplexity,max_iter   the tsne parameters used in visualisation
#' @value none
#' @export


metaClustering<-function(indataframe=NULL,
                         clustername=NULL,
                         metaClustering_method=NULL,
                         usecol=NULL,
                         elbow_test=T,
                         k_value=NULL,
                         #tsne parameters：
                         view_tsne=T,
                         seed=NULL,
                         perplexity=15,
                         max_iter=1500,
                         ...){
  
  if(0){
    
    indataframe=FlowSOM_combined
    usecol=NULL
    #tsne parameters：
    view_tsne=T
    perplexity=15
    max_iter=1500
    
    clustername = "FlowSOM"
    metaClustering_method = "metaClustering_PhenoGraph" #<-- 聚类方法
    k_value=10 #<-- 聚类方法中的k值
    elbow_test=F #<-- 决定是否进行elbow test，当kvalue=NULL时，会自行选择进行elbowtest
    seed=123
  }
        cat("Get expression matrix of cluster centers...\n")
        if (is.null(usecol))
              {usecol<-c(1:ncol(indataframe))}
        clusterparaid<-which(colnames(indataframe)==clustername)
        usecol<-union(usecol, c(clusterparaid))
        cluster_center_expr<-data.frame(indataframe[,usecol]) %>%
              dplyr::group_by_at(clustername) %>%
              dplyr::summarise_if(is.numeric,median)
            
        cluster_abundance<-data.frame(indataframe[,usecol]) %>%
              dplyr::group_by_at(clustername) %>%
              dplyr::summarise(num=n())
            
        cat("MetaClustering using method:",metaClustering_method,"...\n")
        if(elbow_test==T){
               cat("Drawing elow curve...\n")
        }
        #Findvalue<-DetermineNumberOfClusters(data=cluster_center_expr,max=20,kstep = 2,method=metaClustering_method,plot = T)
        
        cc_metacluster<-MetaClustering(data=cluster_center_expr[,-clusterparaid],
                                       method=metaClustering_method,
                                       k_value=k_value,
                                       elbow_test = elbow_test)
        cat("\nMetaclustering cluster centers is finished.\n")
        cat("Start to mapping metaclusters to single cells...\n")
        
        Cluster_arrange<-data.frame(cluster=cluster_center_expr[,clustername],
                                    metacluster=cc_metacluster)
        
        cluster_arrange_fun<-function( cluster_id ){
              cellcluster <- subset(Cluster_arrange,Cluster_arrange[,1]==cluster_id)$metacluster
              return(cellcluster)
            }
        metacluster_result <- apply(as.matrix(indataframe[,colnames(indataframe)==clustername]),1,cluster_arrange_fun)
        
        if(view_tsne==T)
        {
                
                cat("Summarise metacluster information...\n")
                cat("Start to visualise metaclusters with tSNE...\n")
                #tsne (可以调节十多个参数，最重要的两个：perplexity 和 max_iter)
                if(is.null(seed)){
                        seed<-ceiling(runif(1,min=0,max=1)*10000)
                        cat("Seed is not specified, randomly set to: ",seed,"\n")
                        set.seed(seed)
                }else{
                       cat("Seed is set to: ",seed,".\n")
                       set.seed(seed)
                }
          
          
          tsne_result <- Rtsne(cluster_center_expr[,-c(1)], initial_dims = ncol(cluster_center_expr[,-c(1)]),
                               dims = 2, check_duplicates = FALSE, pca = F, perplexity=15,max_iter=1500)$Y
          colnames(tsne_result)<-c("tsne_1","tsne_2")
          
          combine_data_plot<-data.frame(cluster_center_expr,
                                        metacluster=cc_metacluster,
                                        tsne_result,
                                        num=cluster_abundance$num)

          centers<-combine_data_plot %>%
                    dplyr::group_by(metacluster)  %>%
                    dplyr::summarise(tsne_1=median(tsne_1),tsne_2=median(tsne_2))
          
          
          ##做图：visualization using ggplot2
          
          combine_data_plot$metacluster<-as.factor(combine_data_plot$metacluster)
          
          mytheme <- theme(panel.background = element_rect(fill = "white", colour = "black", size = 0.2), #坐标系及坐标轴
                           legend.key = element_rect(fill = "white", colour = "white"), #图标
                           legend.background = (element_rect(colour= "white", fill = "white")))
          
          klab="Cluster number(k_value)"
          if(metaClustering_method=="metaClustering_PhenoGraph") klab="PhenoGraph_k(k_value)"
          ptsnemap<-ggplot(combine_data_plot)+
                    geom_point(aes(x=tsne_1,y=tsne_2,colour=metacluster,size=num),alpha=0.7)+
                    guides(colour = guide_legend(ncol = 2, bycol = T))+
                    scale_size_continuous(range = c(0.1, 5))+
                    labs(title = paste0(klab,": ",k_value))+
                    mytheme+
                    geom_text(data=centers,aes(x=tsne_1,y=tsne_2),label=centers$metacluster,colour="black",size=5)
          
          print(ptsnemap)
        }
      
        cat("Metaclustering is finished successufully.\n")
        return(metacluster_result)
      }
      



#' MetaClustering, elbow_test and tsne visualisation
#'
#' 
#' @param indataframe            event expression matrix or dataframe， containing one clustering channel
#' @param clustername            string,the column name of clustering channel
#' @param metaClustering_method  string,the method of metaclustering to use. could be set to "metaClustering_consensus",
#'                               "metaClustering_hclust","metaClustering_kmeans","metaClustering_PhenoGraph"
#' @param usecol                 numeric, the numeric columns id of indataframe, data in selected columns would be used in metaclustering. 
#' @param elbow_test             TRUE/FALSE, determine whether elbow_test will be carried out before formal metaclustering.
#' @param k_value                numeric, the k value which is used in the metaclustering_method, equal to the number of generated clusters except for PhenoGraph. 
#' @param view_tsne              TRUE/FALSE, determine whether tsne visualisation is performed, use TRUE as default.                               
#' @param seed,perplexity,max_iter   the tsne parameters used in visualisation
#' @value none
#' @export


nodemetaClustering<-function(#indataframe=NULL,
                         #cluster_center_expr=NULL,
                         fSOM=NULL,
                         clustername=NULL,
                         metaClustering_method=NULL,
                         #usecol=NULL,
                         elbow_test=T,
                         k_value=NULL,
                         #tsne parameters：
                         view_tsne=T,
                         seed=NULL,
                         perplexity=15,
                         max_iter=1500,
                         ...){
  
  
  #cat("Get expression matrix of cluster centers...\n")
  cluster_center_expr=fSOM$map$codes
  # if (is.null(usecol))
  # {
  #   #usecol<-c(1:ncol(cluster_center_expr))
  #   usecol<-fSOM$map$colsUsed
  #   }
  # clusterparaid<-which(colnames(indataframe)==clustername)
  # usecol<-union(usecol, c(clusterparaid))
  # cluster_center_expr<-data.frame(indataframe[,usecol]) %>%
  #   dplyr::group_by_at(clustername) %>%
  #   dplyr::summarise_if(is.numeric,median)
  # 
  # cluster_abundance<-data.frame(indataframe[,usecol]) %>%
  #   dplyr::group_by_at(clustername) %>%
  #   dplyr::summarise(num=n())
  
  # cluster_center_expr<-cluster_center_expr[,usecol]
  
  
  cat("MetaClustering using method:",metaClustering_method,"...\n")
  if(elbow_test==T){
    cat("Drawing elow curve...\n")
  }
  #Findvalue<-DetermineNumberOfClusters(data=cluster_center_expr,max=20,kstep = 2,method=metaClustering_method,plot = T)
  
  cc_metacluster<-MetaClustering(data=cluster_center_expr,#[,-clusterparaid],
                                 method=metaClustering_method,
                                 k_value=k_value,
                                 elbow_test = elbow_test)
  cat("\nMetaclustering cluster centers is finished.\n")
 
  
  
  
  
  
  if(view_tsne==T)
  {
    
    cat("Summarise metacluster information...\n")
    cat("Start to visualise metaclusters with tSNE...\n")
    #tsne (可以调节十多个参数，最重要的两个：perplexity 和 max_iter)
    if(is.null(seed)){
      seed<-ceiling(runif(1,min=0,max=1)*10000)
      cat("Seed is not specified, randomly set to: ",seed,"\n")
      set.seed(seed)
    }else{
      cat("Seed is set to: ",seed,".\n")
      set.seed(seed)
    }
    # cluster_center_expr<-fSOM$map$medianValues
    # cluster_center_expr<-cluster_center_expr[ncol(cluster_center_expr) %in% c("Time","Event_length")]
    # cluster_center_expr<-cluster_center_expr[,fSOM$map$colsUsed]
    
    
    tsne_result <- Rtsne(cluster_center_expr[,-c(1)], initial_dims = ncol(cluster_center_expr[,-c(1)]),
                         dims = 2, check_duplicates = FALSE, pca = F, perplexity=15,max_iter=1500)$Y
    colnames(tsne_result)<-c("tsne_1","tsne_2")
    
    combine_data_plot<-data.frame(cluster_center_expr,
                                  metacluster=cc_metacluster,
                                  tsne_result
                                  #num=cluster_abundance$num
                                  )
    
    centers<-combine_data_plot %>%
      group_by(metacluster)  %>%
      summarise(tsne_1=median(tsne_1),tsne_2=median(tsne_2))
    
    
    ##做图：visualization using ggplot2
    
    combine_data_plot$metacluster<-as.factor(combine_data_plot$metacluster)
    
    mytheme <- theme(panel.background = element_rect(fill = "white", colour = "black", size = 0.2), #坐标系及坐标轴
                     legend.key = element_rect(fill = "white", colour = "white"), #图标
                     legend.background = (element_rect(colour= "white", fill = "white")))
    
    klab="Cluster number(k_value)"
    if(metaClustering_method=="metaClustering_PhenoGraph") klab="PhenoGraph_k(k_value)"
    ptsnemap<-ggplot(combine_data_plot)+
      geom_point(aes(x=tsne_1,y=tsne_2,colour=metacluster),alpha=0.7)+
      guides(colour = guide_legend(ncol = 2, bycol = T))+
      #scale_size_continuous(range = c(0.1, 5))+
      labs(title = paste0(klab,": ",k_value))+
      mytheme+
      geom_text(data=centers,aes(x=tsne_1,y=tsne_2),label=centers$metacluster,colour="black",size=5)
    
    print(ptsnemap)
  }
  
  
  
  
  return(cc_metacluster)
}




map_to_singlecells<-function(cc_metacluster,fSOM,view_tsne=T,seed=NULL)

{
  
  cat("Start to mapping metaclusters to single cells...\n")
  
  Cluster_arrange<-data.frame(cluster=c(1:length(cc_metacluster)),
                              metacluster=cc_metacluster
                              )
  

  indataframe<-data.frame(FlowSOM=fSOM$map$mapping[,1])
  
  
  
  cluster_arrange_fun<-function( cluster_id ){
    cellcluster <- subset(Cluster_arrange,Cluster_arrange[,1]==cluster_id)$metacluster
    return(cellcluster)
  }
  metacluster_result <- apply(as.matrix(indataframe[,colnames(indataframe)=="FlowSOM"]),1,cluster_arrange_fun)
  

}
