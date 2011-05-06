#!/usr/bin/env Rscript
### 05.01.11
contigSizes <- read.csv("contig_sizes",header=F,sep="\t")

readData <- function(mappingFile) {
   data <- read.csv(mappingFile,header=F,sep="\t")
   names(data) <- c("chr","start","end","contig","contig_start","contig_end","strand")
   data <- merge(data,contigSizes,by.x="contig",by.y="V1",all.x=T)
   
   data$fraction <- (data$contig_end-data$contig_start)/data$V2
   data$star <- ""
   data[data$fraction >= .5,]$star <- "*"
   data[data$fraction >= .6,]$star <- "**"
   data[data$fraction >= .7,]$star <- "***"
   data[data$fraction >= .8,]$star <- "****"
   data[data$fraction >= .9,]$star <- "*****"
   
   data$text_color <- "black"
   data[data$fraction<.5,]$text_color <- "darkred"
   
   data$color <- "black"
#   data[data$strand=="-",]$color <- "royalblue"
   data
}

blat2dmel <- read.csv("min1000_dmelARMS_3Rinv2contigs.psl",header=F,sep="\t",skip=6,as.is=T)
names(blat2dmel) <- c("match","mismatch","repMatch","Ns","ref_gapCount","ref_gapBases","gapCount","gapBases","strand",
                      "refChr","refSize","refStart","refEnd",
                      "contig","contigLength","start","end",
                      "blocks","blockSizes","ref_blockStarts","blockStarts")


####################################################################################################
pdf("stitchedContigs.pdf",width=5,height=10)
par(mar=c(5,1.5,.5,11),bg="transparent",cex.lab=1,cex.axis=.68,font.main=2,font.lab=2)
for (chrom in c("2L","2R","3L","3Rinverted","4","X")) {

   contigData <- readData(sprintf("chr%s_mapping-1000-500.tsv",chrom))
   contigData <- contigData[order(contigData$start),]
   #row.names(contigData) <- sprintf("%s:%d-%d %s (%d)",contigData$contig,contigData$contig_start,contigData$contig_end,contigData$strand,contigData$start)

   cat(sprintf("Chromosome %s\n",chrom))
   cat(sprintf("Starting with contig hits %d\n",nrow(contigData)))

   contigsSplit <- {}#list(contig=NULL,contigStart=NULL,contigEnd=NULL,refChr=NULL)
   newHits <- {}

   origContigData <- contigData
   for (i in 1:nrow(origContigData)) {
      focal <- origContigData[i,]
      focal$start <- as.numeric(as.vector(focal$start))
      focal$end   <- as.numeric(as.vector(focal$end))
      overlaps <- origContigData[
         (((origContigData$start <= focal$start) & (origContigData$end >= focal$end)) |
          ((origContigData$start >= focal$start) & (origContigData$end <= focal$end)) |
          ((origContigData$start <= focal$start) & (origContigData$end >= focal$start)) |
          ((origContigData$start <= focal$end)   & (origContigData$end >= focal$end))),]
      if (nrow(overlaps)>=5) { 
         cat(sprintf("(%d hits)\t%s %s %d..%d ---",nrow(overlaps),focal$contig,focal$chr,focal$start,focal$end)) 
         #cat(sprintf("(%d hits)\t%s %s %d..%d\n",nrow(overlaps),focal$contig,focal$chr,focal$start,focal$end)) 
         contig2split <- blat2dmel[blat2dmel$contig==focal$contig & blat2dmel$refChr==chrom & blat2dmel$refStart==focal$start & blat2dmel$refEnd==focal$end,]
            
         blockSizes <- as.numeric(unlist(strsplit(contig2split$blockSizes,split=",")))
         ref_blockStarts <- as.numeric(unlist(strsplit(contig2split$ref_blockStarts,split=",")))
         blockStarts <- as.numeric(unlist(strsplit(contig2split$blockStarts,split=",")))
         if (contig2split$strand == "-") { blockStarts <- contig2split$contigLength - blockStarts }

         dist2last <- ref_blockStarts[-1] - ref_blockStarts[-length(ref_blockStarts)]
         where2split <- which(dist2last > 1000)#==max(dist2last))
         where2split <- which(dist2last==max(dist2last))
         cat(sprintf("gap size (%d)",max(dist2last)))
         cat(" gap sizes ")
         cat(dist2last[where2split])
         cat(" where2split indices ")
         cat(where2split)
         cat("\n")

         contigData <- contigData[! (contigData$contig==focal$contig & contigData$start==focal$start & contigData$end==focal$end),]
         newHits2add <- rbind(
            c(chrom,ref_blockStarts[1],ref_blockStarts[where2split]+blockSizes[where2split],
            focal$contig,blockStarts[1],blockStarts[where2split]+blockSizes[where2split],contig2split$strand),
            c(chrom,ref_blockStarts[where2split+1],ref_blockStarts[length(ref_blockStarts)]+blockSizes[length(ref_blockStarts)],
            focal$contig,blockStarts[where2split+1],blockStarts[length(ref_blockStarts)]+blockSizes[length(ref_blockStarts) ],contig2split$strand))

         contigsSplit <- rbind(contigsSplit,list(contig=focal$contig,contigStart=focal$contig_start,contigEnd=focal$contig_end,refChr=focal$chr))
         newHits <- rbind(newHits,
            c(chrom,ref_blockStarts[1],ref_blockStarts[where2split]+blockSizes[where2split],
               as.vector(focal$contig),blockStarts[1],blockStarts[where2split]+blockSizes[where2split],contig2split$strand),
            c(chrom,ref_blockStarts[where2split+1],ref_blockStarts[length(ref_blockStarts)]+blockSizes[length(ref_blockStarts)],
               as.vector(focal$contig),blockStarts[where2split+1],blockStarts[length(ref_blockStarts)]+blockSizes[length(ref_blockStarts) ],contig2split$strand))
      }
#      if (nrow(overlaps) > 1) { contigData[i,]$color <- "red" } 
   }

   newContigData <- contigData
   if (is.null(newHits) == FALSE) {
      newHits <- as.data.frame(newHits)
      names(newHits) <- c("chr","start","end","contig","contig_start","contig_end","strand")
      newHits <- merge(as.data.frame(newHits),contigSizes,by.x="contig",by.y="V1",all.x=T)
      newHits$fraction <- 0
      newHits$star <- "-----"
      newHits$text_color <- "royalblue" 
      newHits$color <- "royalblue" 
      row.names(newHits) <- sprintf("%s:%d-%d %s (%d)",newHits$contig,newHits$contig_start,newHits$contig_end,newHits$strand,newHits$start)

      cat(sprintf("Split apart %d contigs\n",nrow(contigsSplit)))
      cat(sprintf("After  removal %d\n",nrow(contigData)))
      contigData <- rbind(contigData,newHits)
      contigData$start <- as.numeric(contigData$start)
      contigData$end   <- as.numeric(contigData$end)
      cat(sprintf("Final          %d\n",nrow(contigData)))

      newHits2add <- as.data.frame(newHits2add)
      names(newHits2add) <- c("chr","start","end","contig","contig_start","contig_end","strand")
      newContigData <- rbind(contigData[,c("chr","start","end","contig","contig_start","contig_end","strand")],newHits2add)
   }

   contigData <- contigData[order(contigData$start),]
   newContigData <- newContigData[order(newContigData$start),]
   write.table(file=sprintf("chr%s_mapping-1000-500.NEW.tsv",chrom),newContigData,quote=F,sep="\t",row.names=F,col.names=F)

   plot(c(0,0),col="transparent",xlim=c(1,max(contigData$end)),ylim=c(1,nrow(contigData)),axes=F,xlab=chrom,ylab="",yaxs="i")#,xaxs="i")
   axis(1,line=.5)
   mtext(rownames(contigData),side=4,at=1:nrow(contigData),cex=.38,las=2,col=contigData$text_color)
   mtext(contigData$star,side=2,at=1:nrow(contigData),cex=.68,las=2,col=contigData$text_color)
#   points(contigData$start,1:nrow(contigData),pch=15,cex=.5,col=contigData$color,xpd=T)
   segments(contigData$start,1:nrow(contigData),contigData$end,1:nrow(contigData),lwd=2,col=contigData$color,lend=1,xpd=T)
}
dev.off()

####################################################################################################

