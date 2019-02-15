
microscope.counts <- function(input=" ", output=" "){

# Input data
data.in <-  read.csv(input,header=TRUE)                          
data.in <- as.matrix(data.in)

#microscope <- scope
time <- data.in[1,]
time <- as.numeric(time[7:length(time)])
data.in <- data.in[2:nrow(data.in),]

# Initialize data storage
results <- NULL

# Create Output
#create if statement here for the type of experiment (options = bac.counts, pop.counts)

#outfile.id <- paste("./output/",output,".txt", sep="")
#titles.id <- c("time","trt","identifier","org","abundance", "sem")
#write.table(as.matrix(t(titles.id)), file=outfile.id, append=T, row.names=F, col.names=F, sep=",", quote=FALSE)


# Number of fields for different objectives
#NOTE: TURN THIS INTO A IF STATEMENT FOR MICROSCOPE TYPE
x20=4897.85
x40=18783.42
x100=118540.63

identifier <- unique(data.in[,3])

for(i in 1:length(identifier)){
  cur <- data.in[data.in[,3] == identifier[i],]
  
  if(unique(cur[,5])=="Control"){	
    ids <- cur[1,2:6]
    cur <- cur[,7:ncol(cur)]
    fields <- as.numeric(cur[12,])
    fields[fields==20] <- x20
    fields[fields==40] <- x40
    fields[fields==100] <- x100
    cpm <- colMeans(matrix(as.numeric(cur[1:10,]),10,ncol(cur)),na.rm=TRUE)/(as.numeric(cur[11,])*(1-as.numeric(cur[13,])))*1000*fields
    results <- rbind(results,c(ids,cpm))
  }else{
    curS <- cur[cur[,2]=="Syn",]
    idsS <- curS[1,2:6]
    curS <- curS[,7:ncol(curS)]
    fieldsS <- as.numeric(curS[12,])
    fieldsS[fieldsS==20] <- x20
    fieldsS[fieldsS==40] <- x40
    fieldsS[fieldsS==100] <- x100
    cpmS <- colMeans(matrix(as.numeric(curS[1:10,]),10,ncol(curS)),na.rm=TRUE)/(as.numeric(curS[11,])*(1-as.numeric(curS[13,])))*1000*fieldsS
    results <- rbind(results,c(idsS,cpmS))
    
    curP <- cur[cur[,2]=="Phage",]
    idsP <- curP[1,2:6]
    curP <- curP[,7:ncol(curP)]
    fieldsP <- as.numeric(curP[12,])
    fieldsP[fieldsP==20]=x20
    fieldsP[fieldsP==40]=x40
    fieldsP[fieldsP==100]=x100
    cpmP <- colMeans(matrix(as.numeric(curP[1:10,]),10,ncol(curP)),na.rm=TRUE)/(as.numeric(curP[11,])*(1-as.numeric(curP[13,])))*1000*fieldsP
    results <- rbind(results,c(idsP,cpmP))
  }
}

results1 <- results
IDS <- colnames(results[,1:5])
colnames(results) <- c(IDS,time)
outfile.id <- paste("./output/",output,".csv", sep="")
write.csv(results,file=outfile.id,row.names = FALSE)

}