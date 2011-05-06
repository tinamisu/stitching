### 04.19.11
# qsub -l 1hr -cwd -N prepStitching "R CMD BATCH prepStitchCoords.R"

####################################################################################################
readData <- function(datafile) {
   data <- read.csv(datafile,skip=4,as.is=T,sep="\t",header=F)
   # [S1]    [E1]    [S2]    [E2]    [LEN 1] [LEN 2] [% IDY] [LEN R] [LEN Q] [COV R] [COV Q] [TAGS]
   names(data) <- c("start1","end1","start2","end2","len1","len2","percentID","total1","total2","cov1","cov2","name1","name2")
   return(data)
}   

dmel2dsim   <-readData("ref2ref.coords")
dmel2stitch <-readData("stitch2dmel.coords")
dsim2stitch <-readData("stitch2dsim.coords") ### sim 2 stitch
   
contigs <- c("2L","2R","3L","3R","4","X")


### FIGURE 1: stitch vs. (dsim,dmel)
pdf("stitching2ref_dotplot.pdf",height=5,width=5)
par(mar=c(3,4,2,.5),cex.lab=.8,cex.axis=.38,cex.main=1,font.main=2,bg="transparent") # ,mgp=c(1,.5,0)
for (contig in contigs) {
   dmel_data <- dmel2stitch[dmel2stitch$name1==contig & dmel2stitch$name2==contig,]
   dsim_data <- dsim2stitch[dsim2stitch$name1==contig & dsim2stitch$name2==contig,]
   plotRange <- c(1,max(c(dmel_data$total1,dsim_data$total1)))

   plot(c(0,0),col="transparent",xlab="",ylab="",xlim=c(1,max(dmel_data$total2)),ylim=plotRange,main=contig)
   mtext(side=1,line=2,"stitched",col="black")
   mtext(side=2,line=2,"dsim ref",col="blue")
   mtext(side=2,line=3,"dmel ref",col="red")
   lines(plotRange,plotRange,col="gray68")
   segments(dmel_data$start2,dmel_data$start1,dmel_data$end2,dmel_data$end1,lwd=2,lend=1,col="red")#col=gray(1-dmel_data$percentID/100))
   segments(dsim_data$start2,dsim_data$start1,dsim_data$end2,dsim_data$end1,lwd=2,lend=1,col="blue")#col=gray(1-dmel_data$percentID/100))
}
dev.off()



### FIGURE 2: dsim vs. dmel
#./nucmer --prefix=ref_dmel2dsim dmel-all-chromosome-r5.33.fasta dsim-all-chromosome-r1.3.fasta
pdf("dsim2dmel_dotplot.pdf",height=5,width=5)
par(mar=c(3,3,2,.5),cex.lab=.8,cex.axis=.38,cex.main=1,font.main=2,bg="transparent") # ,mgp=c(1,.5,0)
for (contig in contigs) {
   data <- dmel2dsim[dmel2dsim$name1==contig & dmel2dsim$name2==contig,]
   plotRange <- c(1,max(c(data$total1,data$total1)))

   plot(c(0,0),col="transparent",xlab="",ylab="",xlim=c(1,max(data$total2)),ylim=c(1,max(data$total1)),main=contig)
   mtext(side=1,line=2,"dsim ref",col="black")
   mtext(side=2,line=2,"dmel ref",col="black")
   lines(plotRange,plotRange,col="gray68")
   segments(data$start2,data$start1,data$end2,data$end1,lwd=2,lend=1,col="gray28")
}
dev.off()
