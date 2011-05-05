### 04.14.11
### plot

# blatData <- read.csv("contigsVsDsim.psl",header=F,sep="\t",skip=6)
# names(blatData) <- c("match","mismatch","repMatch","Ns","gapCount","gapBases","ref_gapCount","ref_gapBases","strand","contig","contigLength","start","end","refChr","refSize","refStart","refEnd","blocks","blockSizes","blockStarts","ref_blockStarts")
# 
# TEdata <- read.csv("testTE.psl",header=F,sep="\t",skip=6)
# names(TEdata) <- c("match","mismatch","repMatch","Ns","gapCount","gapBases","ref_gapCount","ref_gapBases","strand","contig","contigLength","start","end","refChr","refSize","refStart","refEnd","blocks","blockSizes","blockStarts","ref_blockStarts")
# 
# ### blockSizes,blockStarts,ref_blockStarts
# TEhits <- unlist(lapply(as.vector(TEdata$blockStarts),function(x) { strsplit(x,split=",")}))
# 
# ### order ref chroms by the amount of match to the contig
# 
# ### convert to list
# 
# for (chrom in levels(data$refChr)) {
# 
# yLimits <- {};
# 
# plotData <- data[data$refChr=="3L",]
# plotBlatData <- blatData[blatData$refChr=="3L",]
# 
# yRange <- range(c(unlist(plotData[,c("refStart","refEnd")]),unlist(plotBlatData[,c("refStart","refEnd")])))
# plot(0,0,xlim=c(1,500108),ylim=yRange,col="transparent",xlab="contig",ylab="dsim REF")
# rect(TEdata$start,rep(0,nrow(TEdata)),TEdata$end,rep(yRange[2],nrow(TEdata)),col="yellow",bor="goldenrod")
# arrows(TEdata$start,rep(-5,nrow(TEdata)),TEdata$end,rep(-5,nrow(TEdata)),col="orange")
# segments(plotBlatData$start,plotBlatData$refStart,plotBlatData$end,plotBlatData$refEnd,col="red",lwd=5)
# segments(plotData$start,plotData$refStart,plotData$end,plotData$refEnd)
# 
# #contigName <- "NODE_82870_length_791080_cov_18.549400"


### BAD REGIONS (from MSG)
####################################################################################################
mask <- read.table("../badRegions_msg_S9.1.csv",header=T,as.is=T)


### READ IN BLAT DATA
####################################################################################################
full_blatData <- read.csv("dmelARMS_3Rinv2contigs.psl",header=F,sep="\t",skip=6,as.is=T)
#full_blatData <- read.csv("contigs_63_Dsim.psl",header=F,sep="\t",skip=6,as.is=T)
#blatData <- read.csv("contigs_min790000.psl",header=F,sep="\t",skip=6)
names(full_blatData) <- c("match","mismatch","repMatch","Ns","ref_gapCount","ref_gapBases","gapCount","gapBases","strand","refChr","refSize","refStart","refEnd","contig","contigLength","start","end","blocks","blockSizes","ref_blockStarts","blockStarts")
#names(full_blatData) <- c("match","mismatch","repMatch","Ns","gapCount","gapBases","ref_gapCount","ref_gapBases","strand","contig","contigLength","start","end","refChr","refSize","refStart","refEnd","blocks","blockSizes","blockStarts","ref_blockStarts")

print(paste("starting with:",nrow(full_blatData)))

### FILTER CONTIGS THAT HAVE HITS > 10000
#blatData <- full_blatData[full_blatData$contigLength >= 500000,]
blatData <- full_blatData[full_blatData$contigLength >= 1000,]
#blatData <- full_blatData[full_blatData$contigLength >= 10000,]

#ok <- blatData$match+blatData$mismatch+blatData$repMatch+blatData$Ns < 5000
ok <- blatData$match+blatData$mismatch+blatData$repMatch+blatData$Ns < 500
#ok <- blatData$match+blatData$mismatch+blatData$repMatch+blatData$Ns < 1000
print(paste("removing <1000 bp hits:",length(ok)))
blatData <- blatData[!ok,]
print(paste("retaining:",nrow(blatData)))

### ADD IN TE matches



### SUMMARIZE THE REF CHROMS
####################################################################################################
mainArms <- c("2L","2R","3L","3R","X","4")
plotChroms <- unique(full_blatData[,c("refChr","refSize")])
row.names(plotChroms) <- plotChroms$refChr
plotChroms <- plotChroms[mainArms,]
#plotChroms <- plotChroms[c(mainArms,setdiff(row.names(plotChroms),mainArms)),]
#plotChroms <- plotChroms[order(plotChroms$refSize,decreasing=T),]
plotChroms$plotStarts <- c(0,cumsum(plotChroms$refSize)[-length(plotChroms$refSize)])

#plotContigs <- unique(blatData[,c("contig","contigLength")])
#plotContigs <- plotContigs[order(plotContigs$contigLength,decreasing=T),]
#plotContigs$plotStarts <- c(0,cumsum(plotContigs$contigLength)[-length(plotContigs$contigLength)])

### remove hits from non-arms
blatData <- blatData[blatData$refChr %in% mainArms,]


####################################################################################################
segColors <- c("royalblue","darkred","forestgreen")

pdf("contigs2dsim.pdf",width=20,height=12)
#par(mfrow=c(length(chroms),1),mfrow=c(3,3,0.5,.5))
par(mar=c(14,3,0.5,1),bg="transparent")
yRange <- c(0,sum(plotChroms$refSize))
xRange <- c(0,sum(unique(blatData[,c("contig","contigLength")])$contigLength))

plot(0,0,xlim=xRange,ylim=yRange,xlab="",ylab="",col="transparent",axes=F,xaxs="i",yaxs="i")
#for (m in 1:nrow(mask)) {
#   if (mask[m,]$chr %in% plotChroms$refChr) {
#      rect(xRange[1],mask[m,]$start+plotChroms[plotChroms$refChr==mask[m,]$chr,]$plotStarts,xRange[2],mask[m,]$end+plotChroms[plotChroms$refChr==mask[m,]$chr,]$plotStarts,col="yellow",bor="transparent")
#   }
#}
abline(h=plotChroms$plotStarts,col="gray28")
mtext(plotChroms$refChr,side=2,line=.5,at=plotChroms$plotStarts,las=2)
box(col="gray68")

### order contigs by their position relative to reference
plotContigs <- data.frame(contig=NULL,contigLength=NULL,plotStarts=NULL)
for (refChrom in as.vector(plotChroms$refChr)) {
   contigs2chrom <- blatData[blatData$refChr==refChrom,]
   contigs2chrom <- contigs2chrom[order(contigs2chrom$refStart),]

   if (nrow(plotContigs)==0) { 
      plotContigs <- unique(contigs2chrom[,c("contig","contigLength")])
      plotContigs$plotStarts <- c(0,cumsum(plotContigs$contigLength)[-length(plotContigs$contigLength)])

   } else { 
      addonContigs <- unique(contigs2chrom[contigs2chrom$contig %in% setdiff(contigs2chrom[,"contig"],plotContigs[,"contig"]), c("contig","contigLength")])
      if (nrow(addonContigs) > 0) {
         addonContigs$plotStarts <- c(sum(plotContigs$contigLength),sum(plotContigs$contigLength)+cumsum(addonContigs$contigLength)[-length(addonContigs$contigLength)]) 
         plotContigs <- rbind(plotContigs,addonContigs)
      }
   }

   for (contigName in as.vector(contigs2chrom$contig)) {
      plotData <- blatData[blatData$contig==contigName & blatData$refChr==refChrom,]
      
      ### break up blocks
      blockSizes <- as.numeric(unlist(lapply(as.vector(plotData$blockSizes),function(x) { strsplit(x,split=",")})))
      refStarts <- as.numeric(unlist(lapply(as.vector(plotData$ref_blockStarts),function(x) { strsplit(x,split=",")}))) + plotChroms[plotChroms$refChr==refChrom,]$plotStarts
      
      contigStarts <- as.numeric(unlist(apply(plotData,1,function(y) { 
         x <- as.data.frame(t(y));
         if (x$strand=="-") { as.numeric(as.vector(x$contigLength)) - as.numeric(unlist(strsplit(as.vector(x$blockStarts),split=","))) }
         else { strsplit(as.vector(x$blockStarts),split=",") }
         }))) + plotContigs[plotContigs$contig==contigName,]$plotStarts
      
      contigEnds <- as.numeric(unlist(apply(plotData,1,function(y) { 
         x <- as.data.frame(t(y));
         x_blockSizes <- as.numeric(unlist(strsplit(as.vector(x$blockSizes),split=",")))
         if (x$strand=="-") { as.numeric(as.vector(x$contigLength)) - as.numeric(unlist(strsplit(as.vector(x$blockStarts),split=","))) - x_blockSizes + 1 }
         else { as.numeric(unlist(strsplit(as.vector(x$blockStarts),split=","))) + x_blockSizes }
         }))) + plotContigs[plotContigs$contig==contigName,]$plotStarts
      
      segments(contigStarts,refStarts,contigEnds,refStarts+blockSizes,col="royalblue",lend=1,lwd=2)
   }

}

#abline(v=plotContigs$plotStarts,col="gray68")
#mtext(plotContigs$contig,side=1,line=.5,at=plotContigs$plotStarts,las=2,cex=.68)

dev.off()
#save.image()
