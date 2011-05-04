### 04.19.11
# qsub -l 1hr -cwd -N prepStitching "R CMD BATCH prepStitchCoords.R"

plotBlocks <- function(plotData,yPos,plotColor) {
   blockSizes <- as.numeric(unlist(lapply(as.vector(plotData$blockSizes),function(x) { strsplit(x,split=",")})))

   contigStarts <- as.numeric(unlist(apply(plotData,1,function(y) { 
         x <- as.data.frame(t(y));
         if (x$strand=="-") { as.numeric(as.vector(x$contigLength)) - as.numeric(unlist(strsplit(as.vector(x$blockStarts),split=","))) }
         else { strsplit(as.vector(x$blockStarts),split=",") }
         })))
   contigEnds <- as.numeric(unlist(apply(plotData,1,function(y) { 
         x <- as.data.frame(t(y));
         x_blockSizes <- as.numeric(unlist(strsplit(as.vector(x$blockSizes),split=",")))
         if (x$strand=="-") { as.numeric(as.vector(x$contigLength)) - as.numeric(unlist(strsplit(as.vector(x$blockStarts),split=","))) - x_blockSizes + 1 }
         else { as.numeric(unlist(strsplit(as.vector(x$blockStarts),split=","))) + x_blockSizes }
         })))

   segments(contigStarts,rep(yPos,length(blockSizes)),contigEnds,rep(yPos,length(blockSizes)),col=plotColor,lend=1,lwd=3)
}


      
#   refStarts <- as.numeric(unlist(lapply(as.vector(plotData$ref_blockStarts),function(x) { strsplit(x,split=",")}))) #+ plotChroms[plotChroms$refChr==refChrom,]$plotStarts
#   contigStarts <- as.numeric(unlist(apply(plotData,1,function(y) { 
#         x <- as.data.frame(t(y));
#         if (x$strand=="-") { as.numeric(as.vector(x$contigLength)) - as.numeric(unlist(strsplit(as.vector(x$blockStarts),split=","))) }
#         else { strsplit(as.vector(x$blockStarts),split=",") }
#         }))) #+ plotContigs[plotContigs$contig==contigName,]$plotStarts
#      
#   contigEnds <- as.numeric(unlist(apply(plotData,1,function(y) { 
#         x <- as.data.frame(t(y));
#         x_blockSizes <- as.numeric(unlist(strsplit(as.vector(x$blockSizes),split=",")))
#         if (x$strand=="-") { as.numeric(as.vector(x$contigLength)) - as.numeric(unlist(strsplit(as.vector(x$blockStarts),split=","))) - x_blockSizes + 1 }
#         else { as.numeric(unlist(strsplit(as.vector(x$blockStarts),split=","))) + x_blockSizes }
#         }))) #+ plotContigs[plotContigs$contig==contigName,]$plotStarts
#   segments(contigStarts,refStarts,contigEnds,refStarts+blockSizes,col=plotColor,lend=1,lwd=2)


min_match <- 500
min_contig_size <- 1000

### READ IN BLAT DATA
####################################################################################################
blat2dmel <- read.csv("contigs_63_Dmel.psl",header=F,sep="\t",skip=6,as.is=T)
names(blat2dmel) <- c("match","mismatch","repMatch","Ns","gapCount","gapBases","ref_gapCount","ref_gapBases","strand","contig","contigLength","start","end","refChr","refSize","refStart","refEnd","blocks","blockSizes","blockStarts","ref_blockStarts")
blat2dmel$alnLength <- blat2dmel$match + blat2dmel$mismatch + blat2dmel$repMatch + blat2dmel$Ns

blat2dsim <- read.csv("contigs_63_Dsim.psl",header=F,sep="\t",skip=6,as.is=T)
names(blat2dsim) <- c("match","mismatch","repMatch","Ns","gapCount","gapBases","ref_gapCount","ref_gapBases","strand","contig","contigLength","start","end","refChr","refSize","refStart","refEnd","blocks","blockSizes","blockStarts","ref_blockStarts")
blat2dsim$alnLength <- blat2dsim$match + blat2dsim$mismatch + blat2dsim$repMatch + blat2dsim$Ns

contigs <- read.csv("contig_sizes",sep="\t",as.is=T,header=F)
names(contigs) <- c("contig","size")

mapped_contigs <- read.csv("all_mapping-1000-500.tsv",header=F,sep="\t",as.is=T)
names(mapped_contigs) <- c("chr","chrStart","chrEnd","contig","start","end","strand")

### FILTER DATA
####################################################################################################
# [1] "Dmel starting with: 202117"
cat(sprintf("Starting with %d contigs",nrow(contigs)))
contigs <- contigs[contigs$size>=min_contig_size,]
cat(sprintf("%d contigs > %d",nrow(contigs),min_contig_size))
cat(sprintf("%d contigs > %d found in dmel",length(intersect(contigs$contig,blat2dmel$contig)),min_contig_size))
cat(sprintf("%d contigs > %d found in dsim",length(intersect(contigs$contig,blat2dsim$contig)),min_contig_size))
cat(sprintf("%d contigs > %d found in dmel & dsim",length(intersect(contigs$contig,intersect(blat2dmel$contig,blat2dsim$contig))),min_contig_size))

pdf("contig2dmeldsim.pdf",height=3,width=5)
par(mar=c(3,3.8,2,3),cex.lab=.8,cex.axis=.38,cex.main=.68,mgp=c(1,.1,0))
for (c in 1:nrow(contigs)) {
   contig <- contigs[c,]

   contig <- contigs[contigs$contig=="NODE_40389_length_138103_cov_18.340912",]

   data_dmel <- blat2dmel[blat2dmel$contig==contig$contig & blat2dmel$alnLength>=min_match,]
   data_dsim <- blat2dsim[blat2dsim$contig==contig$contig & blat2dsim$alnLength>=min_match,]

   data_dmel <- data_dmel[order(data_dmel$alnLength),]
   data_dsim <- data_dsim[order(data_dsim$alnLength),]

   hits2dmel <- nrow(data_dmel)
   hits2dsim <- nrow(data_dsim)
   contig_title <- 'unmapped'
   if (contig$contig %in% mapped_contigs$contig) {
      contig_title <- paste(mapped_contigs[mapped_contigs$contig==contig$contig,c("chr","chrStart","chrEnd")],collapse=" ",sep="\n")
      #contig_title <- sprintf("%s:%s..%s",as.vector(unlist(mapped_contigs[mapped_contigs$contig==contig$contig,c("chr","chrStart","chrEnd")])))

      plot(

   }

   plot(c(1,contigs[c,]$size),c(0,0),type="l",lwd=1,lend=1,xlim=c(1,contig$size),ylim=c(0,hits2dmel+hits2dsim),
        main=contig_title,xlab=paste(contig$contig,contig$size),ylab="",axes=F)
   abline(h=hits2dmel+.5,col="gray68",lty=2)
   axis(1)

   #rect(focal_dmel$start,1,focal_dmel$end,hits2dmel+hits2dsim+1,col="gray88",bor="transparent")
         
   if (hits2dmel>0) {
      segments(data_dmel$start,1:hits2dmel,data_dmel$end,1:hits2dmel,lwd=1,lend=1)#,col=1:hits2dmel)
      #segments(data_dmel$start,1:hits2dmel,data_dmel$end,1:hits2dmel,lwd=1,lend=1,col=1:hits2dmel)
      for (j in 1:hits2dmel) { plotBlocks(data_dmel[j,],j,j) }
      mtext(sprintf("%s:%d",data_dmel$refChr,data_dmel$refStart),side=2,line=.05,at=c(1:hits2dmel),cex=.28,las=2)
   }
   
   if (hits2dsim>0) {
      segments(data_dsim$start,c(1:hits2dsim)+hits2dmel,data_dsim$end,c(1:hits2dsim)+hits2dmel,lwd=1,lend=1)#,col=1:hits2dsim)
      for (j in 1:hits2dsim) { plotBlocks(data_dsim[j,],j+hits2dmel,j) }
      mtext(sprintf("%s:%d",data_dsim$refChr,data_dsim$refStart),side=4,line=.05,at=c(1:hits2dsim)+hits2dmel,cex=.28,las=2)
   }
}
dev.off()
