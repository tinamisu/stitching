#!/usr/bin/env Rscript
source("~/R/getopts.R")

### 04.14.11
### plot BLAT hits (after filtering)
opts <- getopts()
blatFile <- opts$f
species <- opts$s

### BAD REGIONS (from MSG)
####################################################################################################
#mask <- read.table("../badRegions_msg_S9.1.csv",header=T,as.is=T)
#segColors <- c("royalblue","darkred","forestgreen")

### READ IN BLAT DATA
####################################################################################################
full_blatData <- read.csv(blatFile,header=F,sep="\t",skip=6,as.is=T)
names(full_blatData) <- c("match","mismatch","repMatch","Ns","ref_gapCount","ref_gapBases","gapCount","gapBases","strand","refChr","refSize","refStart","refEnd","contig","contigLength","start","end","blocks","blockSizes","ref_blockStarts","blockStarts")
cat(sprintf("starting with: %d\n",nrow(full_blatData)))

### FILTER CONTIGS THAT HAVE HITS > 10000
blatData <- full_blatData[full_blatData$contigLength >= 1000,]
ok <- blatData$match+blatData$mismatch+blatData$repMatch+blatData$Ns < 500
cat(sprintf("removing <1000 bp hits: %d\n",length(ok)))
blatData <- blatData[!ok,]
cat(sprintf("retaining: %d\n",nrow(blatData)))


### SUMMARIZE THE REF CHROMS
####################################################################################################
mainArms <- c("2L","2R","3L","3R","3Rinverted","X","4")
plotChroms <- unique(full_blatData[,c("refChr","refSize")])
row.names(plotChroms) <- plotChroms$refChr
useableChr <- intersect(mainArms,as.vector(plotChroms$refChr))

plotChroms <- plotChroms[useableChr,]
#plotChroms <- plotChroms[c(mainArms,setdiff(row.names(plotChroms),mainArms)),]
#plotChroms <- plotChroms[order(plotChroms$refSize,decreasing=T),]
plotChroms$plotStarts <- c(0,cumsum(plotChroms$refSize)[-length(plotChroms$refSize)])

### remove hits from non-arms
blatData <- blatData[blatData$refChr %in% useableChr,]
cat("ref chroms mapped to: ")
cat(useableChr)
cat("\n")


####################################################################################################
bitmap(file=sprintf("contigs2%s.bmp",species),bg="transparent", width=100, height=100)
par(mar=c(1,3,0.5,1),bg="transparent")
yRange <- c(0,sum(plotChroms$refSize))
xRange <- c(0,sum(unique(blatData[,c("contig","contigLength")])$contigLength))

plot(0,0,xlim=xRange,ylim=yRange,xlab="",ylab="",col="transparent",axes=F,xaxs="i",yaxs="i")
abline(h=plotChroms$plotStarts,col="gray28")
mtext(plotChroms$refChr,side=2,line=.5,at=plotChroms$plotStarts,las=2)
box(col="gray68")

### order contigs by their position relative to reference
plotContigs <- data.frame(contig=NULL,contigLength=NULL,plotStarts=NULL)
for (refChrom in useableChr) {
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
      refStarts <- as.numeric(unlist(apply(plotData,1,function(y) { 
         x <- as.data.frame(t(y));
         if (x$strand=="-") { as.numeric(as.vector(x$refSize)) - as.numeric(unlist(strsplit(as.vector(x$ref_blockStarts),split=","))) }
         else { strsplit(as.vector(x$ref_blockStarts),split=",") }
         }))) + plotChroms[plotChroms$refChr==refChrom,]$plotStarts
      
      refEnds <- as.numeric(unlist(apply(plotData,1,function(y) { 
         x <- as.data.frame(t(y));
         x_blockSizes <- as.numeric(unlist(strsplit(as.vector(x$blockSizes),split=",")))
         if (x$strand=="-") { as.numeric(as.vector(x$refSize)) - as.numeric(unlist(strsplit(as.vector(x$ref_blockStarts),split=","))) - x_blockSizes + 1 }
         else { as.numeric(unlist(strsplit(as.vector(x$ref_blockStarts),split=","))) + x_blockSizes }
         }))) + plotChroms[plotChroms$refChr==refChrom,]$plotStarts

      contigStarts <- as.numeric(unlist(lapply(as.vector(plotData$blockStarts),function(x) { strsplit(x,split=",")}))) + plotContigs[plotContigs$contig==contigName,]$plotStarts
      
      segments(contigStarts,refStarts,contigStarts+blockSizes,refEnds,col="royalblue",lend=1,lwd=2)
   }

}
#abline(v=plotContigs$plotStarts,col="gray68")
#mtext(plotContigs$contig,side=1,line=.5,at=plotContigs$plotStarts,las=2,cex=.68)
rug(plotContigs$plotStarts,ticksize=-0.01,side=1,col="gray28")

dev.off()
#save.image()
