#!/usr/bin/env Rscript

### 04.14.11
### filter and plot BLAT hits
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
cat(sprintf("starting with: %d",nrow(full_blatData)))

### FILTER CONTIGS THAT HAVE HITS > 10000
blatData <- full_blatData[full_blatData$contigLength >= 1000,]
ok <- blatData$match+blatData$mismatch+blatData$repMatch+blatData$Ns < 500
cat(sprintf("removing <1000 bp hits: %d",length(ok)))
blatData <- blatData[!ok,]
cat(sprintf("retaining: %d",nrow(blatData)))


### SUMMARIZE THE REF CHROMS
####################################################################################################
mainArms <- c("2L","2R","3L","3R","X","4")
plotChroms <- unique(full_blatData[,c("refChr","refSize")])
row.names(plotChroms) <- plotChroms$refChr
plotChroms <- plotChroms[mainArms,]
#plotChroms <- plotChroms[c(mainArms,setdiff(row.names(plotChroms),mainArms)),]
#plotChroms <- plotChroms[order(plotChroms$refSize,decreasing=T),]
plotChroms$plotStarts <- c(0,cumsum(plotChroms$refSize)[-length(plotChroms$refSize)])

### remove hits from non-arms
blatData <- blatData[blatData$refChr %in% as.vector(plotChroms$refChr),]


####################################################################################################
pdf(sprintf("contigs2%s.pdf",species),width=20,height=12)
par(mar=c(14,3,0.5,1),bg="transparent")
yRange <- c(0,sum(plotChroms$refSize))
xRange <- c(0,sum(unique(blatData[,c("contig","contigLength")])$contigLength))

plot(0,0,xlim=xRange,ylim=yRange,xlab="",ylab="",col="transparent",axes=F,xaxs="i",yaxs="i")
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
