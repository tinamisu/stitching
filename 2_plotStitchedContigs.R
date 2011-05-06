#!/usr/bin/env Rscript
### 05.01.11
contigSizes <- read.csv("contig_sizes",header=F,sep="\t")
chromSizes  <- read.csv("dmel_size.tsv",header=F,sep="\t")
row.names(chromSizes) <- as.vector(chromSizes$V1)

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
#   data[data$fraction<.5,]$text_color <- "darkblue"
   
   data$color <- "black"
#   data[data$strand=="-",]$color <- "royalblue"
   data
}


####################################################################################################
for (chrom in c("2L","2R","3L","3Rinverted","4","X")) {
   contigData <- readData(sprintf("chr%s_mapping-1000-500.tsv",chrom))
   contigData <- contigData[order(contigData$start),]

   row.names(contigData) <- sprintf("%s %d-%d (%d)",contigData$contig,contigData$contig_start,contigData$contig_end,contigData$start)
   #row.names(contigData) <- sprintf("%s:%d-%d %s (%d)",contigData$contig,contigData$contig_start,contigData$contig_end,contigData$strand,contigData$start)

   tmp <- as.data.frame(cbind(contigData$contig,contigData$contig_start,contigData$contig_end))
   changeColors <- which(duplicated(tmp) | duplicated(tmp,fromLast = TRUE))
   if (length(changeColors)>0) { contigData[changeColors,]$text_color <- "royalblue" }

   cat(sprintf("Chromosome %s\n",chrom))
   cat(sprintf("Starting with contig hits %d\n",nrow(contigData)))

   contigData <- contigData[order(contigData$start),]

   ####################################################################################################
   pdf(sprintf("stitchedContigs_chr%s.pdf",chrom),width=5,height=18)
   par(mar=c(5,1.5,.5,11),bg="transparent",cex.lab=1,cex.axis=.68,font.main=2,font.lab=2)
   plot(c(0,0),col="transparent",xlim=c(1,max(contigData$end)),ylim=c(1,nrow(contigData)),axes=F,xlab=chrom,ylab="",yaxs="i")#,xaxs="i")
   axis(1,line=.5)
   mtext(rownames(contigData),side=4,at=1:nrow(contigData),cex=.28,las=2,col=contigData$text_color)
   mtext(contigData$star,side=2,at=1:nrow(contigData),cex=.68,las=2,col=contigData$text_color)
#   points(contigData$start,1:nrow(contigData),pch=15,cex=.5,col=contigData$color,xpd=T)
   segments(contigData$start,1:nrow(contigData),contigData$end,1:nrow(contigData),lwd=2,col=contigData$color,lend=1,xpd=T)
   dev.off()

   ####################################################################################################
   contig2start <- unique(contigData[,c("contig","V2")])
   contig2start$plotStart <- c(0,cumsum(contig2start$V2)[-nrow(contig2start)])

   bitmap(file=sprintf("stitchedContigsMap%s.bmp",chrom),bg="transparent", width=150, height=100)
   par(mar=c(1,5,.5,.5),bg="transparent")
   xRange <- c(0,chromSizes[chrom,]$V2)
   yRange <- c(0,max(contig2start$plotStart)+contig2start[nrow(contig2start),]$V2)

   plot(0,0,xlim=xRange,ylim=yRange,xlab="",ylab="",col="transparent",axes=F,xaxs="i",yaxs="i")
   box(col="gray28")
   abline(h=contig2start$plotStart,col="gray68")
   mtext(contig2start$contig,side=2,line=.1,at=contig2start$plotStart,las=2,cex=.28)
   mtext(sprintf("dmel %s",chrom),side=1,line=.005,font=2,cex=.8)

   for (c in 1:nrow(contig2start)) {
      plotData <- contigData[contigData$contig==as.vector(contig2start[c,]$contig),]
      segments(plotData$start,plotData$contig_start+contig2start[c,]$plotStart,plotData$end,plotData$contig_end+contig2start[c,]$plotStart,
         col="red",lend=1,lwd=3)
   }
   dev.off()
}


### DETERMINE CONTIGS THAT MAPPED FULLY
####################################################################################################
system("cat chr*_mapping-1000-500.tsv > all_mapping-1000-500.tsv")
mapped <- readData("all_mapping-1000-500.tsv")

contigCounts <- as.data.frame(table(mapped$contig))
chimeras <- contigCounts[contigCounts$Freq>1,]

mapped$missing <- mapped$V2-(mapped$contig_end-mapped$contig_start)
mapped[mapped$contig %in% chimeras$Var1,]$missing <- NA
mapped$missing_fraction <- mapped$missing / mapped$V2
#mapped[is.na(mapped$missing_fraction)==F & mapped$missing_fraction>.1,]

pdf("uniqueMapping_fraction.pdf",width=5,height=3)
par(mar=c(5,5,1,.5),bg="transparent")
hist(mapped$missing_fraction)
plot(mapped$V2,mapped$missing,type="p",col="gray68")
dev.off()

cat(sprintf("Total Number of contigs mapped uniquely %d\n",nrow(contigCounts)))
cat(sprintf("Number of contigs involved in chimeras %d\n",nrow(chimeras)))
write.table(mapped[is.na(mapped$missing_fraction)==F & mapped$missing_fraction <= .1,],file="unique_contigs_full.tsv",sep="\t",row.names=F,col.names=F,quote=F)

cat(sprintf("Number possible chimera but not mapped %d\n",nrow(mapped[is.na(mapped$missing_fraction)==F & mapped$missing_fraction > .1,])))
write.table(mapped[is.na(mapped$missing_fraction)==F & mapped$missing_fraction > .1,],file="possible_chimera_contigs.tsv",sep="\t",row.names=F,col.names=F,quote=F)
write.table(chimeras,file="chimera_contigs.tsv",sep="\t",row.names=F,col.names=F,quote=F)
