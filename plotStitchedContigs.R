data <- read.csv("all_mapping-1000-500.tsv",header=F,sep="\t")
names(data) <- c("chr","start","end","contig","contig_start","contig_end","strand")
contigSizes <- read.csv("contig_sizes",header=F,sep="\t")
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

data$color <- "red"
data[data$strand=="-",]$color <- "blue"

pdf("stitchedContigs.pdf",width=8,height=28)
par(mar=c(5,1.5,.5,11),bg="transparent",cex.lab=1,cex.axis=.68,font.main=2,font.lab=2)
for (contig in c("2L","2R","3L","3R","4","X")) {
   contigData <- data[data$chr==contig,]
   contigData <- contigData[order(contigData$start),]
   row.names(contigData) <- sprintf("%s:%d-%d %s",contigData$contig,contigData$contig_start,contigData$contig_end,contigData$strand)
   plot(c(0,0),col="transparent",xlim=c(1,max(contigData$end)),ylim=c(1,nrow(contigData)),axes=F,xlab=contig,ylab="",yaxs="i")#,xaxs="i")
   axis(1,line=.5)
   mtext(rownames(contigData),side=4,at=1:nrow(contigData),cex=.38,las=2,col=contigData$text_color)
   mtext(contigData$star,side=2,at=1:nrow(contigData),cex=.68,las=2,col=contigData$text_color)
   segments(contigData$start,1:nrow(contigData),contigData$end,1:nrow(contigData),lwd=2,col=contigData$color,lend=1,xpd=T)
   #segments(contigData$start,contigData$contig_start,contigData$end,contigData$contig_end,lwd=2,col=contigData$color)
}
dev.off()
