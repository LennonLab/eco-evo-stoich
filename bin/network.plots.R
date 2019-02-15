inf.network <- function(type = "", cid = "", 
                        phage.color="", bac.color="", 
                        phage.lty = "", bac.lty = ""){
  require(igraph)
  if(type == "cstat"){
    #Read in data and modify form
    data.in <- read.csv(paste("../output/infmatrix",cid,".csv",sep = ""), header = F)
    
    data.in[1,1] <- 999
    data.in[,1] <- as.numeric(as.character(data.in[,1]))
    data.in <- as.matrix(data.in)
    
  } else {
    data.in <- read.csv(paste("../output/infmatrix",cid,".csv",sep = ""), header = F)
    
    data.in[1,1] <- 999
    data.in[,1] <- as.numeric(as.character(data.in[,1]))
    data.in <- as.matrix(data.in)   
  }


  rownames(data.in)<-c()
  colnames(data.in)<-c()
  
  # grab out time vectors
  bacteria.time <- data.in[1,2:(dim(data.in))[[2]]]
  phage.time <- as.numeric(data.in[2:dim(data.in)[[1]],1])
  #phage.time <- c(0,9,23,72,129,166)
  #bacteria.time <- c(-6,9,23,72,100,129,166)
  
  #Trim data
  data.in <- data.in[2:dim(data.in)[[1]], 2:(dim(data.in))[[2]]]
  
  #generate a square adjacency matrix, creating a place to store the data
  mat <- matrix(0,nrow=sum(dim(data.in)),ncol=sum(dim(data.in)))
  
  nrows <- dim(data.in)[[1]]
  big.nrows <- dim(mat)[[2]]
  
  # create the phage-bacteria edges between time points
  mat[1:nrows,I(nrows+1):big.nrows] <- data.in
  
  # convert adjacency matrix into graph object
  g <- graph.adjacency(mat, mode = "directed", weighted = T)
  
  # extract weights
  #wts <- g[[9]][[4]]$weight
  wts <- E(g)$weight
  
  # create layout matrix
  locs <- cbind(c(phage.time,bacteria.time),
                c(rep(8,length(phage.time)),rep(1,length(bacteria.time))))
  
  # to account for rescaling of time into -1 to 1, we need to set up axes labels appropriately
  xs <- c(0,40,80,120,160)
  scaled.xs <- 2*(xs - min(locs[,1]))/(max(locs[,1]) - min(locs[,1])) - 1
  scaled.infect <- 2*(locs[,1] - min(locs[,1]))/(max(locs[,1]) - min(locs[,1])) - 1
  
  # modify vertex labels
  vlabels <- c(phage.time,bacteria.time)
  
  
  #tiff(filename = paste("./figures/inf.ts",cid,".tif",sep = ""),
  #width = 480, height = 640, units = px, pointsize = 12,
  #bg = "white", family = "sans", type = "cairo")
  
  plot.inf <- plot(g,vertex.shape = "circle",
                      vertex.label = NA,
                      #vertex.label = vlabels,
                      #vertex.label.color = c(rep("white",length(phage.time)),
                      #                       rep(bac.color,length(bacteria.time))),
                      #vertex.label.font = 4, vertex.label.cex = 1,
                      #vertex.size = c(pvsize*60,bvsize*20),
                      vertex.color = c(rep(phage.color,length(phage.time)),
                                       rep("white",length(bacteria.time))),
                      edge.color = "black", edge.width = 8* wts, 
                      edge.arrow.mode = 0,layout = locs,
                      #xlab = list("Time (d)",cex = 1.25),
                   asp = 0,ylim = c(-1,1),
                      main = cid)
  
  # draw scaled axes
  axis(1, at = scaled.xs,labels = xs,cex.axis = 1,font = 1)
  axis(2, at = c(-1,1),labels = c("host","phage"),cex.axis = 1.25,font = 1, las = 1)
}


inf.network.ts <- function(type = "", cid = "", 
                           phage.color="", bac.color="",
                           phage.lty = "", bac.lty = ""){
  require(igraph)
  require(plotrix)
  if(type == "cstat"){
    #Read in data and modify form
    data.in <- read.csv(paste("../output/infmatrix",cid,".csv",sep = ""), header = F)
    
    data.in[1,1] <- 999
    data.in[,1] <- as.numeric(as.character(data.in[,1]))
    data.in <- as.matrix(data.in)
    
    tsdat <- read.csv("../output/cid-means.csv")
    tsdat <- subset(tsdat, tsdat$day > -11 & tsdat$cID == cid)
    tsdat <- na.omit(tsdat)
    
    # rescale time series variables
    dat.y.B <- log10(tsdat$abd[tsdat$microbe == "Syn"])
    dat.x.B <- tsdat$day[tsdat$microbe == "Syn"]
    #dat.x.B <- dat.x.B[-which(is.na(dat.y.B))]
    #dat.y.B <- dat.y.B[-which(is.na(dat.y.B))]
    
    dat.y.P <- log10(tsdat$abd[tsdat$microbe == "Phage"])
    dat.x.P <- tsdat$day[tsdat$microbe == "Phage"]
    #dat.x.P <- dat.x.P[-which(is.na(dat.y.P))]
    #dat.y.P <- dat.y.P[-which(is.na(dat.y.P))]
    
  } else {
    data.in <- read.csv(paste("../output/infmatrix",cid,".csv",sep = ""), header = F)
    
    data.in[1,1] <- 999
    data.in[,1] <- as.numeric(as.character(data.in[,1]))
    data.in <- as.matrix(data.in)
    
    # need to import file with same format as previous data file
    tsdat <- read.csv("../output/cstat-means.csv")
    tsdat <- subset(tsdat, tsdat$day > -11)
    
    if(cid == "N"){
      # rescale time series variables
      dat.y.B <- log10(tsdat$NmeanSI)
      dat.x.B <- tsdat$day
      #dat.x.B <- dat.x.B[-which(is.na(dat.y.B))]
      #dat.y.B <- dat.y.B[-which(is.na(dat.y.B))]
      
      dat.y.P <- log10(tsdat$NmeanP)
      dat.x.P <- tsdat$day
      dat.x.P <- dat.x.P[-which(is.na(dat.y.P))]
      dat.y.P <- dat.y.P[-which(is.na(dat.y.P))]
      
      dat.P <- cbind(dat.x.P,dat.y.P)
      dat.P <- na.omit(dat.P)
      
      } else {
      dat.y.B <- log10(tsdat$PmeanSI)
      dat.x.B <- tsdat$day
        #dat.x.B <- dat.x.B[-which(is.na(dat.y.B))]
        #dat.y.B <- dat.y.B[-which(is.na(dat.y.B))]
        
      dat.y.P <- log10(tsdat$PmeanP)
      dat.x.P <- tsdat$day
      dat.x.P <- dat.x.P[-which(is.na(dat.y.P))]
      dat.y.P <- dat.y.P[-which(is.na(dat.y.P))]      
    }
    

  }

  
  rownames(data.in)<-c()
  colnames(data.in)<-c()
  
  # grab out time vectors
  bacteria.time <- data.in[1,2:(dim(data.in))[[2]]]
  phage.time <- as.numeric(data.in[2:dim(data.in)[[1]],1])
  #phage.time <- c(0,9,23,72,129,166)
  #bacteria.time <- c(-6,9,23,72,100,129,166)
  
  #Trim data
  data.in <- data.in[2:dim(data.in)[[1]], 2:(dim(data.in))[[2]]]
  
  #generate a square adjacency matrix, creating a place to store the data
  mat <- matrix(0,nrow=sum(dim(data.in)),ncol=sum(dim(data.in)))
  
  nrows <- dim(data.in)[[1]]
  big.nrows <- dim(mat)[[2]]
  
  # create the phage-bacteria edges between time points
  mat[1:nrows,I(nrows+1):big.nrows] <- data.in
  
  # convert adjacency matrix into graph object
  g <- graph.adjacency(mat, mode = "directed", weighted = T)
  
  # extract weights
  #wts <- g[[9]][[4]]$weight
  wts <- E(g)$weight
  
  # create layout matrix
  locs <- cbind(c(phage.time,bacteria.time),
              c(rep(8,length(phage.time)),rep(1,length(bacteria.time))))
  
  # to account for rescaling of time into -1 to 1, we need to set up axes labels appropriately
  xs <- c(0,40,80,120,160)
  scaled.xs <- 2*(xs - min(locs[,1]))/(max(locs[,1]) - min(locs[,1])) - 1
  scaled.infect <- 2*(locs[,1] - min(locs[,1]))/(max(locs[,1]) - min(locs[,1])) - 1
  
  # modify vertex labelsS
  vlabels <- c(phage.time,bacteria.time)
  
  
  if(cid == "N"){
    plot.inf.ts <- plot(g,vertex.shape = "circle",
                      vertex.label = NA,
                      #vertex.label = vlabels,
                      #vertex.label.color = c(rep("white",length(phage.time)),
                      #                       rep(bac.color,length(bacteria.time))),
                      #vertex.label.font = 4, vertex.label.cex = 1,
                      vertex.size = 10,
                      vertex.color = c(rep(phage.color,length(phage.time)),
                                        rep("white",length(bacteria.time))),
                      edge.color = "black", edge.width = 6* wts, 
                      edge.arrow.mode = 0,layout = locs,
                      xlab = list("Time (d)",cex = 1.5),asp = 0,ylim = c(-10,1),
                      main = list("N-limited", cex = 1.5, font = 2))
  } else {
      if (cid == "P"){
        plot.inf.ts <- plot(g,vertex.shape = "circle",
                        vertex.label = NA,
                        #vertex.label = vlabels,
                        #vertex.label.color = c(rep("white",length(phage.time)),
                        #                       rep(bac.color,length(bacteria.time))),
                        #vertex.label.font = 4, vertex.label.cex = 1,
                        vertex.size = 10,
                        vertex.color = c(rep(phage.color,length(phage.time)),
                                         rep("white",length(bacteria.time))),
                        edge.color = "black", edge.width = 6* wts, 
                        edge.arrow.mode = 0,layout = locs,
                        xlab = list("Time (d)",cex = 1.5),asp = 0,ylim = c(-10,1),
                        main = list("P-limited", cex = 1.5, font = 2))   
      } else {
        plot.inf.ts <- plot(g,vertex.shape = "circle",
                        vertex.label = NA,
                        #vertex.label = vlabels,
                        #vertex.label.color = c(rep("white",length(phage.time)),
                        #                       rep(bac.color,length(bacteria.time))),
                        #vertex.label.font = 4, vertex.label.cex = 1,
                        vertex.size = 10,
                        vertex.color = c(rep(phage.color,length(phage.time)),
                                         rep("white",length(bacteria.time))),
                        edge.color = "black", edge.width = 6* wts, 
                        edge.arrow.mode = 0,layout = locs,
                        xlab = list("Time (d)",cex = 1.5),asp = 0,ylim = c(-10,1),
                        main = list(cid, cex = 1.5, font = 2))
  }
}
  
  # draw scaled axes
  axis(1, at = scaled.xs,labels = xs,cex.axis = 1.25,font = 1)
  #axis(2, at = c(-1,1),labels = c("cyanobacteria","phage"),cex.axis = 1,font = 1)
  
  
  dat.y.B.scaled <- 6*(dat.y.B-min(dat.y.B))/(max(dat.y.B)-min(dat.y.B))-10
  dat.y.P.scaled <- 6*(dat.y.P-min(dat.y.P))/(max(dat.y.P)-min(dat.y.P))-8.4
  dat.x.B.scaled <- 2*(dat.x.B-min(locs[,1]))/(max(locs[,1])-min(locs[,1]))-1
  dat.x.P.scaled <- 2*(dat.x.P-min(locs[,1]))/(max(locs[,1])-min(locs[,1]))-1
  
  # plot time series
  lines(dat.y.B.scaled~dat.x.B.scaled,lwd=3,col=bac.color,lty = bac.lty)
  lines(dat.y.P.scaled~dat.x.P.scaled,lwd=3,col=phage.color, lty = phage.lty)
  
  
  ##for ts data and evolution plot
  ticks <- c(5,6,7,8,9)
  labels <- sapply(ticks, function(i) as.expression(bquote(10^ .(i))))
  axis(2,at=c(-10,-8,-6,-4,-2,-1,1),labels=c(labels,"host","phage"), las = 1 , cex.axis = 1.25)
  #axis.break(2, -1.5, style = 'slash')
  
  
  # add reference lines
  i <- 1
  
  for(i in 1:length(bacteria.time)){
    xloc <- scaled.infect[I(length(phage.time)+i)]
    yloc <- dat.y.B.scaled[dat.x.B.scaled == xloc]
    
    yloc <- min(c(dat.y.B.scaled[dat.x.B.scaled == xloc],
                  dat.y.P.scaled[dat.x.P.scaled == xloc]))
    
    segments(xloc,yloc,xloc,-1.3,lty=2)	
  }
  
  # Show start date of phage addition
  #segments(dat.x.scaled[dat.x==0],dat.y.scaled[dat.x==0]-0.5,
           #dat.x.scaled[dat.x==0],dat.y.scaled[dat.x==0]+0.5,lwd=3)
}

