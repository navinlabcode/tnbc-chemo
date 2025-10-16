ElbowPoint_gives_optK <- function(k, v) {
  df  <- data.frame(k = k, v = v)
  i <- pathviewr::find_curve_elbow(df)
  return(k[i])
}

PAC_gives_optK <- function(results) {
  maxK <- length(results)
  Kvec <- 2:maxK
  x1 <- 0.1
  x2 <- 0.9 # threshold defining the intermediate sub-interval
  PAC <- rep(NA, length(Kvec))
  names(PAC) <- paste("K=", Kvec, sep = "") # from 2 to maxK
  for (i in Kvec) {
    M <- results[[i]]$consensusMatrix
    Fn <- ecdf(M[lower.tri(M)])
    PAC[i - 1] <- Fn(x2) - Fn(x1)
  } # end for i
  # The optimal K
  optK <- Kvec[which.min(PAC)]
  return(optK)
}

ConsensusClusterPlus <- function( d=NULL,
                                  maxK = 3,
                                  reps=10,
                                  pItem=0.8,
                                  pFeature=1,
                                  clusterAlg="hc",
                                  title="untitled_consensus_cluster",
                                  innerLinkage="average",
                                  finalLinkage="average",
                                  distance="pearson",
                                  ml=NULL,
                                  tmyPal=NULL,
                                  seed=NULL,
                                  plot=NULL,
                                  writeTable=FALSE,
                                  weightsItem=NULL,
                                  weightsFeature=NULL,
                                  verbose=F,
				  corUse="everything" ) {

  ##description: runs consensus subsamples 
  if(is.null(seed)==TRUE){
    seed=timeSeed = as.numeric(Sys.time())
  }
  set.seed(seed)

  #distance=ifelse( inherits(d,"dist"), attr( d, "method" ), "pearson" )


  if(is.null(ml)==TRUE){
    
    if ( ! any(class( d ) %in% c( "dist", "matrix", "ExpressionSet" )) ) {
      stop("d must be a matrix, distance object or ExpressionSet (eset object)")
    }

    if ( inherits( d, "dist" ) ) {
      ## if d is a distance matrix, fix a few things so that they don't cause problems with the analysis
      ##  Note, assumption is that if d is a distance matrix, the user doesn't want to sample over the row features
      if ( is.null( attr( d, "method" ) ) ) {
        attr( d, "method" ) <- distance <- "unknown - user-specified"
      }
      if ( is.null( distance ) || ( distance != attr( d, "method" ) ) ) {
        distance <- attr( d, "method" ) 
      }
      
      if ( ( ! is.null( pFeature ) ) && ( pFeature < 1 ) ) {
        message( "Cannot use the pFeatures parameter when specifying a distance matrix as the data object, setting pFeature to 1.\n" )
        pFeature <- 1
      }
      if ( ! is.null( weightsFeature ) ) {
        message( "Cannot use the weightsFeature parameter when specifying a distance matrix as the data object\n" )
        weightsFeature <- NULL
      }

    } else {
      if ( is.null( distance ) ) {
        ## we should never get here, but just in case
	message("no specified distance, setting distance to Pearson");
        distance <- "pearson"
      }
    }

    if ( ( clusterAlg == "km" ) && inherits( distance, "character" ) && ( distance != "euclidean" ) ) {
      message( "Note: The km (kmeans) option only supports a euclidean distance metric when supplying a data matrix.  If you want to cluster a distance matrix, use a different algorithm such as 'hc' or 'pam'.  Changing distance to euclidean")
      distance <- 'euclidean'
    }


    if ( inherits( d,"ExpressionSet" ) ) {
      d <- exprs(d)
    }

    ml <- ccRun( d=d,
                 maxK=maxK,
                 repCount=reps,
                 diss=inherits(d,"dist"),
                 pItem=pItem,
                 pFeature=pFeature,
                 innerLinkage=innerLinkage,
                 clusterAlg=clusterAlg,
                 weightsFeature=weightsFeature,
                 weightsItem=weightsItem,
                 distance=distance,
                 verbose=verbose,
		 corUse=corUse)
  }
  res=list();
  
  ##make results directory
  if((is.null(plot)==FALSE | writeTable) & !file.exists(paste(title,sep=""))){
    dir.create(paste(title,sep=""))
  }
  
  ##write log file
  log <- matrix( ncol=2,
                 byrow=T,
                 c("title",title,
                   "maxK",maxK,
                   "input matrix rows",ifelse ( inherits( d, "matrix" ), nrow(d), "dist-mat" ), 
                   "input matrix columns",ifelse ( inherits( d, "matrix" ), ncol(d), ncol( as.matrix(d) ) ), 
                   "number of bootstraps",reps,
                   "item subsampling proportion",pItem,
                   "feature subsampling proportion",ifelse( is.null(pFeature), 1, pFeature ),
                   "cluster algorithm",clusterAlg,
                   "inner linkage type",innerLinkage,
                   "final linkage type",finalLinkage,
                   "correlation method",distance,
                   "plot",if(is.null(plot)) NA else plot,
                   "seed",if(is.null(seed)) NA else seed))
  colnames(log) = c("argument","value")
  if(writeTable){
    write.csv(file=paste(title,"/",title,".log.csv",sep=""), log,row.names=F)
  }
  if(is.null(plot)){
    ##nothing
  }else if(plot=="pngBMP"){
    bitmap(paste(title,"/","consensus%03d.png",sep=""))
  }else if(plot=="png"){
    png(paste(title,"/","consensus%03d.png",sep=""))

  }else if (plot=="pdf"){
    pdf(onefile=TRUE, paste(title,"/","consensus.pdf",sep=""))
  }else if (plot=="ps"){
    postscript(onefile=TRUE, paste(title,"/","consensus.ps",sep=""))
  }	
  
  colorList=list()
  colorM = rbind() #matrix of colors.
  
  #18 colors for marking different clusters
  thisPal <- c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F","#FF7F00","#CAB2D6","#6A3D9A","#FFFF99","#B15928",
               "#bd18ea", #magenta
               "#2ef4ca", #aqua
               "#f4cced", #pink,
               "#f4cc03", #lightorange
               "#05188a", #navy,
               "#e5a25a", #light brown
               "#06f106", #bright green
               "#85848f", #med gray
               "#000000", #black
               "#076f25", #dark green
               "#93cd7f",#lime green
               "#4d0776", #dark purple
               "#ffffff" #white
               )
  
  ##plot scale
  colBreaks=NA
  if(is.null(tmyPal)==TRUE){
    colBreaks=10
    tmyPal = myPal(colBreaks)
  }else{
    colBreaks=length(tmyPal)
  }
  sc = cbind(seq(0,1,by=1/( colBreaks) )); rownames(sc) = sc[,1]
  sc = cbind(sc,sc)
  heatmap(sc, Colv=NA, Rowv=NA, symm=FALSE, scale='none', col=tmyPal, na.rm=TRUE,labRow=rownames(sc),labCol=F,main="consensus matrix legend")

  for (tk in 2:maxK){
    if(verbose){
      message(paste("consensus ",tk))
    }
    fm = ml[[tk]]
    hc=hclust( as.dist( 1 - fm ), method=finalLinkage);
    message("clustered")	
    ct = cutree(hc,tk)
    names(ct) = colnames(d)
    if(any(class(d)=="dist")){
      names(ct) = colnames(as.matrix(d))
    }
    c = fm

    colorList = setClusterColors(res[[tk-1]][[3]],ct,thisPal,colorList)
    pc = c
    pc=pc[hc$order,] #pc is matrix for plotting, same as c but is row-ordered and has names and extra row of zeros.

    

    if(!is.null(plot) && plot=="pngBMP"){
         pc = pc[,hc$order ] #mod for no tree
        pc = rbind(pc,0)
	#no dendrogram if pngBMP
	oc = colorList[[1]][hc$order] #mod for no tree
	heatmap(pc, Colv = NA, Rowv = NA, symm = FALSE, scale = "none", col = tmyPal, na.rm = TRUE, labRow = F, labCol = F, mar = c(5, 5), main = paste("consensus matrix k=", 
                tk, sep = ""), ColSideCol = oc)
    }else{
        pc = rbind(pc,0)
    #former with tree:
    	heatmap(pc, Colv=as.dendrogram(hc), Rowv=NA, symm=FALSE, scale='none', col=tmyPal, na.rm=TRUE,labRow=F,labCol=F,mar=c(5,5),main=paste("consensus matrix k=",tk,sep="") , ColSideCol=colorList[[1]])
    }

    legend("topright",legend=unique(ct),fill=unique(colorList[[1]]),horiz=FALSE )

    res[[tk]] = list(consensusMatrix=c,consensusTree=hc,consensusClass=ct,ml=ml[[tk]],clrs=colorList)
    colorM = rbind(colorM,colorList[[1]]) 
  }
  CDF_curv_values <- CDF(ml)
  clusterTrackingPlot(colorM[,res[[length(res)]]$consensusTree$order])
  if(is.null(plot)==FALSE){
    dev.off();
  }
  saveRDS(CDF_curv_values, file = paste0(title,"/","CDF_values.rds")) ### YY: Export the CDF values to identify the elbow point ### 
  res[[1]] = colorM
  if(writeTable){
    for(i in 2:length(res)){
      write.csv(file=paste(title,"/",title,".k=",i,".consensusMatrix.csv",sep=""), res[[i]]$consensusMatrix)
      write.table(file=paste(title,"/",title,".k=",i,".consensusClass.csv",sep=""), res[[i]]$consensusClass,col.names = F,sep=",")
    }
  }
  return(res)
}


calcICL = function(res,title="untitled_consensus_cluster",plot=NULL,writeTable=FALSE){
  #calculates and plots cluster consensus and item consensus
  cc=rbind()
  cci = rbind()
  sumRes=list()
  colorsArr=c()
  
  #make results directory
  if((is.null(plot)==FALSE | writeTable) & !file.exists(paste(title,sep=""))){
	dir.create(paste(title,sep=""))
  }
  if(is.null(plot)){
    #to screen
  }else if(plot=="pdf"){
    pdf(onefile=TRUE, paste(title,"/","icl.pdf",sep=""))
  }else if(plot=="ps"){
    postscript(onefile=TRUE, paste(title,"/","icl.ps",sep=""))
  }else if (plot=="png"){
    png(paste(title,"/","icl%03d.png",sep=""))
  }else if (plot=="pngBMP"){
    bitmap(paste(title,"/","icl%03d.png",sep=""))
  }

  par(mfrow=c(3,1),mar=c(4,3,2,0))

  for (k in 2:length(res)){ #each k
    eiCols = c();
    o = res[[k]]
    m = o$consensusMatrix
    m = triangle(m,mode=2)
    for (ci in sort(unique(o$consensusClass))){ #each cluster in k
	items = which(o$consensusClass==ci)
	nk = length(items)
	mk = sum( m[items,items], na.rm=T)/((nk*(nk-1))/2)
	cc=rbind(cc,c(k,ci,mk)) #cluster-consensus
	
      for (ei in rev(res[[2]]$consensusTree$order) ){
		denom = if (ei %in% items) { nk - 1} else { nk }
        	mei = sum( c(m[ei,items],m[items,ei]), na.rm=T)/denom  # mean item consensus to a cluster.
		cci = rbind(cci,c(k,ci,ei,mei)) #cluster, cluster index, item index, item-consensus
      }
      eiCols = c(eiCols, rep(ci,length(o$consensusClass)) )
    }
	  
	  cck = cci[which(cci[,1]==k),] #only plot the new k data.

	  #group by item, order by cluster i
	  w=lapply(split(cck,cck[,3]), function(x) { y=matrix(unlist(x),ncol=4); y[order(y[,2]),4] }) 
	  q = matrix(as.numeric(unlist(w)),ncol=length(w),byrow=F)
	  q = q[,res[[2]]$consensusTree$order] #order by leave order of k=2
 	  #q is a matrix of k rows and sample columns, values are item consensus of sample to the cluster.

	  thisColors = unique(cbind(res[[k]]$consensusClass,res[[k]]$clrs[[1]]))
	  thisColors=thisColors[order(as.numeric(thisColors[,1])),2]
	  colorsArr=c(colorsArr,thisColors)
	  sumRes[[k]] = rankedBarPlot(q,thisColors,cc=res[[k]]$consensusClass[res[[2]]$consensusTree$order],paste("k=",k,sep="") )
  }

  ys=cs=lab=c()
  lastk=cc[1,1]
  for(i in 1:length(colorsArr)){
    if(lastk != cc[i,1]){
      ys=c(ys,0,0)
      cs=c(cs,NA,NA)
      lastk=cc[i,1]
      lab=c(lab,NA,NA)
    }
    ys=c(ys,cc[i,3])
    cs=c(cs,colorsArr[i])
    lab=c(lab,cc[i,1])
  }
  names(ys) = lab
  par(mfrow=c(3,1),mar=c(4,3,2,0))
  barplot(ys,col=cs,border=cs,main="cluster-consensus",ylim=c(0,1),las=1)
  if(is.null(plot)==FALSE){
	  dev.off()
  }
  colnames(cc) = c("k","cluster","clusterConsensus")
  colnames(cci) = c("k","cluster","item","itemConsensus")
  cci[,"item"] = names(res[[2]]$consensusClass)[ cci[,"item"] ]
  #type cci
  cci = data.frame( k=as.numeric(cci[,"k"]), cluster=as.numeric(cci[,"cluster"]), item=cci[,"item"], itemConsensus=as.numeric(cci[,"itemConsensus"])) 
  
  #write to file.
  if(writeTable){
	write.csv(file=paste(title,"/",title,".summary.cluster.consensus.csv",sep=""),row.names=F, cc)
	write.csv(file=paste(title,"/",title,".summary.item.consensus.csv",sep=""), row.names=F, cc)
  }
  return(list(clusterConsensus=cc,itemConsensus=cci))
}


ccRun <- function( d=d,
                   maxK=NULL,
                   repCount=NULL,
                   diss=inherits( d, "dist" ),
                   pItem=NULL,
                   pFeature=NULL,
                   innerLinkage=NULL,
                   distance=NULL, #ifelse( inherits(d,"dist"), attr( d, "method" ), "euclidean" ),@@@@@
                   clusterAlg=NULL,
                   weightsItem=NULL,
                   weightsFeature=NULL,
                   verbose=NULL,
		   corUse=NULL) {
  m = vector(mode='list', repCount)
  ml = vector(mode="list",maxK)
  n <- ifelse( diss, ncol( as.matrix(d) ), ncol(d) )
  mCount = mConsist = matrix(c(0),ncol=n,nrow=n)
  ml[[1]] = c(0);

  if (is.null( distance ) ) distance <- 'euclidean'  ## necessary if d is a dist object and attr( d, "method" ) == NULL

  acceptable.distance <- c( "euclidean", "maximum", "manhattan", "canberra", "binary","minkowski",
                      "pearson", "spearman" )
  
  main.dist.obj <- NULL
  if ( diss ){
    main.dist.obj <- d
    ## reset the pFeature & weightsFeature params if they've been set (irrelevant if d is a dist matrix)
    if ( ( !is.null(pFeature) ) &&
         ( pFeature < 1 ) ) {
      message( "user-supplied data is a distance matrix; ignoring user-specified pFeature parameter\n" )
      pFeature <- 1 # set it to 1 to avoid problems with sampleCols
    }
    if ( ! is.null( weightsFeature ) ) {
      message( "user-supplied data is a distance matrix; ignoring user-specified weightsFeature parameter\n" )
      weightsFeature <- NULL  # set it to NULL to avoid problems with sampleCols
    }
  } else { ## d is a data matrix
    ## we're not sampling over the features
    if ( ( clusterAlg != "km" ) &&
         ( is.null( pFeature ) ||
           ( ( pFeature == 1 ) && is.null( weightsFeature ) ) ) ) {
      ## only generate a main.dist.object IFF 1) d is a matrix, 2) we're not sampling the features, and 3) the algorithm isn't 'km'
      if ( inherits( distance, "character" ) ) {
        if ( ! distance %in%  acceptable.distance  &  ( class(try(get(distance),silent=T))!="function") ) stop("unsupported distance.")

	if(distance=="pearson" | distance=="spearman"){
          main.dist.obj <- as.dist( 1-cor(d,method=distance,use=corUse ))
        }else if( class(try(get(distance),silent=T))=="function"){
          main.dist.obj <- get(distance)( t( d )   )
	}else{
          main.dist.obj <- dist( t(d), method=distance )
        }
        attr( main.dist.obj, "method" ) <- distance  
      } else stop("unsupported distance specified.")
    } else {
      ## pFeature < 1 or a weightsFeature != NULL
      ## since d is a data matrix, the user wants to sample over the gene features, so main.dist.obj is left as NULL
    }
  }
 

  for (i in 1:repCount){
    if(verbose){
      message(paste("random subsample",i));
    }
    ## take expression matrix sample, samples and genes
    sample_x = sampleCols( d, pItem, pFeature, weightsItem, weightsFeature )

    this_dist = NA
    if ( ! is.null( main.dist.obj ) ) {
      boot.cols <- sample_x$subcols
      this_dist <- as.matrix( main.dist.obj )[ boot.cols, boot.cols ]
      if ( clusterAlg != "km" ) {
        ## if this isn't kmeans, then convert to a distance object
        this_dist <- as.dist( this_dist )
        attr( this_dist, "method" ) <- attr( main.dist.obj, "method" )
      }
    } else {
      ## if main.dist.obj is NULL, then d is a data matrix, and either:
      ##   1) clusterAlg is 'km'
      ##   2) pFeatures < 1 or weightsFeatures have been specified, or
      ##   3) both
      ## so we can't use a main distance object and for every iteration, we will have to re-calculate either
      ##   1) the distance matrix (because we're also sampling the features as well), or
      ##   2) the submat (if using km) 

      if ( clusterAlg != "km" )  {
        if ( ! distance %in% acceptable.distance &  ( class(try(get(distance),silent=T))!="function")  ) stop("unsupported distance.")
	if( ( class(try(get(distance),silent=T))=="function") ){
          this_dist <- get(distance)( t( sample_x$submat ) )
	}else{
	  if( distance == "pearson" | distance == "spearman"){
            this_dist <- as.dist( 1-cor(sample_x$submat,use=corUse,method=distance) )
	  }else{
            this_dist <- dist( t( sample_x$submat ), method= distance  )
	  }
	}
        attr( this_dist, "method" ) <- distance  
      } else {
        ## if we're not sampling the features, then grab the colslice
        if ( is.null( pFeature ) ||
            ( ( pFeature == 1 ) && is.null( weightsFeature ) ) ) {
          this_dist <- d[, sample_x$subcols ]
        } else {
          if ( is.na( sample_x$submat ) ) {
            stop( "error submat is NA" )
          }
          
          this_dist <- sample_x$submat
        } 
      }
    }
                  
    ## cluster samples for HC.
    this_cluster=NA
    if(clusterAlg=="hc"){
      this_cluster = hclust( this_dist, method=innerLinkage)
    }
    ##mCount is possible number of times that two sample occur in same random sample, independent of k
    ##mCount stores number of times a sample pair was sampled together.
    mCount <- connectivityMatrix( rep( 1,length(sample_x[[3]])),
                                  mCount,
                                  sample_x[[3]] ) 

    ##use samples for each k		
    for (k in 2:maxK){
      if(verbose){
        message(paste("  k =",k))
      }
      if (i==1){
        ml[[k]] = mConsist #initialize
      }
      this_assignment=NA
      if(clusterAlg=="hc"){
        ##prune to k for hc
        this_assignment = cutree(this_cluster,k)

      
      }else if(clusterAlg=="km"){
        ##this_dist should now be a matrix corresponding to the result from sampleCols
        this_assignment <- kmeans( t( this_dist ),
                                   k,
                                   iter.max = 10,
                                   nstart = 1,
                                   algorithm = c("Hartigan-Wong") )$cluster
      }else if ( clusterAlg == "pam" ) {
        this_assignment <- pam( x=this_dist,
                                k,
                                diss=TRUE,
                                metric=distance, 
                                cluster.only=TRUE )
      } else{
        ##optional cluterArg Hook.
        this_assignment <- get(clusterAlg)(this_dist, k)
      }
      ##add to tally				
      ml[[k]] <- connectivityMatrix( this_assignment,
                                     ml[[k]],
                                     sample_x[[3]] )
    }
  }
	

  ##consensus fraction
  res = vector(mode="list",maxK)
  for (k in 2:maxK){
    ##fill in other half of matrix for tally and count.
    tmp = triangle(ml[[k]],mode=3)
    tmpCount = triangle(mCount,mode=3)
    res[[k]] = tmp / tmpCount
    res[[k]][which(tmpCount==0)] = 0
  }
  message("end fraction")
  return(res)
}


connectivityMatrix <- function( clusterAssignments, m, sampleKey){
  ##input: named vector of cluster assignments, matrix to add connectivities
  ##output: connectivity matrix
  names( clusterAssignments ) <- sampleKey 
  cls <- lapply( unique( clusterAssignments ), function(i) as.numeric( names( clusterAssignments[ clusterAssignments %in% i ] ) ) )  #list samples by clusterId

  for ( i in 1:length( cls ) ) {
    nelts <- 1:ncol( m )
    cl <- as.numeric( nelts %in% cls[[i]] ) ## produces a binary vector
    updt <- outer( cl, cl ) #product of arrays with * function; with above indicator (1/0) statement updates all cells to indicate the sample pair was observed int the same cluster;
    m <- m + updt
  }
  return(m)
}



sampleCols <- function( d,
                        pSamp=NULL,
                        pRow=NULL,
                        weightsItem=NULL,
                        weightsFeature=NULL ){
## returns a list with the sample columns, as well as the sub-matrix & sample features (if necessary)
##  if no sampling over the features is performed, the submatrix & sample features are returned as NAs
##  to reduce memory overhead


  space <- ifelse( inherits( d, "dist" ), ncol( as.matrix(d) ), ncol(d) )
  sampleN <- floor(space*pSamp)
  sampCols <- sort( sample(space, sampleN, replace = FALSE, prob = weightsItem) )

  this_sample <- sampRows <- NA
  if ( inherits( d, "matrix" ) ) {
    if ( (! is.null( pRow ) ) &&
         ( (pRow < 1 ) || (! is.null( weightsFeature ) ) ) ) {
      ## only sample the rows and generate a sub-matrix if we're sampling over the row/gene/features
      space = nrow(d)
      sampleN = floor(space*pRow)
      sampRows = sort( sample(space, sampleN, replace = FALSE, prob = weightsFeature) )
      this_sample <- d[sampRows,sampCols]
      dimnames(this_sample) <- NULL
    } else {
      ## do nothing
    }
  }
  return( list( submat=this_sample,
                subrows=sampRows,
                subcols=sampCols ) )
}

CDF=function(ml,breaks=100){
  #plot CDF distribution
  plot(c(0),xlim=c(0,1),ylim=c(0,1),col="white",bg="white",xlab="consensus index",ylab="CDF",main="consensus CDF", las=2)
  k=length(ml)
  this_colors = rainbow(k-1)
  areaK = c()
  for (i in 2:length(ml)){
    v=triangle(ml[[i]],mode=1)

    #empirical CDF distribution. default number of breaks is 100    
    h = hist(v, plot=FALSE, breaks=seq(0,1,by=1/breaks))
    h$counts = cumsum(h$counts)/sum(h$counts)

    #calculate area under CDF curve, by histogram method.
    thisArea=0
    for (bi in 1:(length(h$breaks)-1)){
       thisArea = thisArea + h$counts[bi]*(h$breaks[bi+1]-h$breaks[bi]) #increment by height by width
       bi = bi + 1
    }
    areaK = c(areaK,thisArea)
    lines(h$mids,h$counts,col=this_colors[i-1],lwd=2,type='l')
  }
  legend(0.8,0.5,legend=paste(rep("",k-1),seq(2,k,by=1),sep=""),fill=this_colors)

  #plot area under CDF change.
  deltaK=areaK[1] #initial auc at k=2
  for(i in 2:(length(areaK))){
    #proportional increase relative to prior K.
    deltaK = c(deltaK,( areaK[i] - areaK[i-1])/areaK[i-1])
  }
  plot(1+(1:length(deltaK)),y=deltaK,xlab="k",ylab="relative change in area under CDF curve",main="Delta area",type="b")
  return(invisible(deltaK))
}

CDF_give_deltaK=function(ml,breaks=100){
  k=length(ml)
  this_colors = rainbow(k-1)
  areaK = c()
  for (i in 2:length(ml)){
    v=triangle(ml[[i]]['consensusMatrix'],mode=1)

    #empirical CDF distribution. default number of breaks is 100    
    h = hist(v, plot=FALSE, breaks=seq(0,1,by=1/breaks))
    h$counts = cumsum(h$counts)/sum(h$counts)

    #calculate area under CDF curve, by histogram method.
    thisArea=0
    for (bi in 1:(length(h$breaks)-1)){
       thisArea = thisArea + h$counts[bi]*(h$breaks[bi+1]-h$breaks[bi]) #increment by height by width
       bi = bi + 1
    }
    areaK = c(areaK,thisArea)
  }

  #plot area under CDF change.
  deltaK=areaK[1] #initial auc at k=2
  for(i in 2:(length(areaK))){
    #proportional increase relative to prior K.
    deltaK = c(deltaK,( areaK[i] - areaK[i-1])/areaK[i-1])
  }
  return(deltaK)
}


myPal = function(n=10){
  #returns n colors
  seq = rev(seq(0,255,by=255/(n)))
  palRGB = cbind(seq,seq,255)
  rgb(palRGB,maxColorValue=255)
}

setClusterColors = function(past_ct,ct,colorU,colorList){
	#description: sets common color of clusters between different K
	newColors = c()
	if(length(colorList)==0){
		#k==2
		newColors = colorU[ct]
		colori=2
	}else{
		newColors = rep(NULL,length(ct))
		colori = colorList[[2]]
		mo=table(past_ct,ct)
		m=mo/apply(mo,1,sum)
			for(tci in 1:ncol(m)){ # for each cluster
				maxC = max(m[,tci])
				pci = which(m[,tci] == maxC)				
				if( sum(m[,tci]==maxC)==1 & max(m[pci,])==maxC & sum(m[pci,]==maxC)==1  )  {
				#if new column maximum is unique, same cell is row maximum and is also unique
				##Note: the greatest of the prior clusters' members are the greatest in a current cluster's members.
					newColors[which(ct==tci)] = unique(colorList[[1]][which(past_ct==pci)]) # one value
				}else{ #add new color.
					colori=colori+1
					newColors[which(ct==tci)] = colorU[colori]
				}
			}
	}
	return(list(newColors,colori,unique(newColors) ))
}

clusterTrackingPlot = function(m){
  #description: plots cluster tracking plot
  #input: m - matrix where rows are k, columns are samples, and values are cluster assignments.
  plot(NULL,xlim=c(-0.1,1),ylim=c(0,1),axes=FALSE,xlab="samples",ylab="k",main="tracking plot")
  for(i in 1:nrow(m)){
    rect(  xleft=seq(0,1-1/ncol(m),by=1/ncol(m)),  ybottom=rep(1-i/nrow(m),ncol(m)) , xright=seq(1/ncol(m),1,by=1/ncol(m)), ytop=rep(1-(i-1)/nrow(m),ncol(m)), col=m[i,],border=NA)   
  }
  #hatch lines to indicate samples
  xl = seq(0,1-1/ncol(m),by=1/ncol(m))
  segments(  xl, rep(-0.1,ncol(m)) , xl, rep(0,ncol(m)), col="black")    #** alt white and black color?
  ypos = seq(1,0,by=-1/nrow(m))-1/(2*nrow(m))
  text(x=-0.1,y=ypos[-length(ypos)],labels=seq(2,nrow(m)+1,by=1))
}

triangle = function(m,mode=1){
  #mode=1 for CDF, vector of lower triangle.
  #mode==3 for full matrix.
  #mode==2 for calcICL; nonredundant half matrix coun
  #mode!=1 for summary 
  n=dim(m)[1]
  nm = matrix(0,ncol=n,nrow=n)
  fm = m


  nm[upper.tri(nm)] = m[upper.tri(m)] #only upper half
  
  fm = t(nm)+nm
  diag(fm) = diag(m)
  
  nm=fm
  nm[upper.tri(nm)] = NA
  diag(nm) = NA
  vm = m[lower.tri(nm)]
  
  if(mode==1){
    return(vm) #vector 		
  }else if(mode==3){
    return(fm) #return full matrix
  }else if(mode == 2){
    return(nm) #returns lower triangle and no diagonal. no double counts.
  }
  
}


rankedBarPlot=function(d,myc,cc,title){
	colors = rbind() #each row is a barplot series
	byRank = cbind()

	spaceh = 0.1 #space between bars
	for(i in 1:ncol(d)){
	  byRank = cbind(byRank,sort(d[,i],na.last=F))
	  colors = rbind(colors,order(d[,i],na.last=F))
	}
	maxH = max(c(1.5,apply(byRank,2,sum)),na.rm=T) #maximum height of graph
	
	#barplot largest to smallest so that smallest is in front.
	barp = barplot( apply(byRank,2,sum) ,  col=myc[colors[,1]] ,space=spaceh,ylim=c(0,maxH),main=paste("item-consensus", title),border=NA,las=1  )
	for(i in 2:nrow(byRank)){
	  barplot( apply(matrix(byRank[i:nrow(byRank),],ncol=ncol(byRank))  ,2,sum), space=spaceh,col=myc[colors[,i]],ylim=c(0,maxH), add=T,border=NA,las=1  )
	}
	xr=seq(spaceh,ncol(d)+ncol(d)*spaceh,(ncol(d)+ncol(d)*spaceh)/ncol(d)  )
	#class labels as asterisks
	text("*",x=xr+0.5,y=maxH,col=myc[cc],cex=1.4) #rect(xr,1.4,xr+1,1.5,col=myc[cc] )
}
