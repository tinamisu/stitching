### 04.29.11
contig_sizes <- read.csv("contig_sizes",header=F,sep="\t")

mapped <- read.csv("all_mapping-1000-500.tsv",header=F,sep="\t",as.is=T)
names(mapped) <- c("refChr","refStart","refEnd","contig","start","end","strand")
contigCounts <- as.data.frame(table(mapped$contig))
chimeras <- contigCounts[contigCounts$Freq>1,]

#mapped$contigLength <- as.numeric(unlist(lapply(strsplit(mapped$contig,split="_"),"[[",4)))
mapped_info <- merge(mapped,contig_sizes,by.x="contig",by.y="V1",all.x=T)
mapped_info$missing <- mapped_info$V2-(mapped_info$end-mapped_info$start)
mapped_info[mapped_info$contig %in% chimeras$Var1,]$missing <- NA
mapped_info$missing_fraction <- mapped_info$missing / mapped_info$V2
mapped_info[is.na(mapped_info$missing_fraction)==F & mapped_info$missing_fraction>.1,]

pdf("uniqueMapping_fraction.pdf",width=5,height=3)
par(mar=c(5,5,1,.5),bg="transparent")
hist(mapped_info$missing_fraction)
plot(mapped_info$V2,mapped_info$missing,type="p",col="gray68")
dev.off()

cat(sprintf("Total Number of contigs mapped uniquely %d",nrow(contigCounts)))
cat(sprintf("Number of contigs involved in chimeras %d",nrow(chimeras)))
write.table(mapped_info[is.na(mapped_info$missing_fraction)==F & mapped_info$missing_fraction <= .1,],file="unique_contigs_full.tsv",sep="\t",row.names=F,col.names=F,quote=F)

cat(sprintf("Number possible chimera but not mapped %d",nrow(mapped_info[is.na(mapped_info$missing_fraction)==F & mapped_info$missing_fraction > .1,])))
write.table(mapped_info[is.na(mapped_info$missing_fraction)==F & mapped_info$missing_fraction > .1,],file="possible_chimera_contigs.tsv",sep="\t",row.names=F,col.names=F,quote=F)
write.table(chimeras,file="chimera_contigs.tsv",sep="\t",row.names=F,col.names=F,quote=F)
                                                             
#print(mapped_info[mapped_info$missing_fraction>.9 & is.na(mapped_info$missing_fraction)==F,])
#blat2dmel[blat2dmel$contig=="NODE_42593_length_9788_cov_18.991112",]



### Everything hit between c(3874907,17560827) on dmel needs to be reversed and reordered!
inv3R_bounds <- c(3874901,17560829)
invInternal <- mapped[mapped$refChr=="3R" & (mapped$refEnd > inv3R_bounds[1] & mapped$refStart < inv3R_bounds[2]),]
#invExternal <- mapped[mapped$refChr=="3R" & (mapped$refEnd < 3874901 | mapped$refStart > 17560829),]

invInternal$newRefStart <- inv3R_bounds[1] + (inv3R_bounds[2]-invInternal$refStart+1)
invInternal$newRefEnd   <- inv3R_bounds[1] + (inv3R_bounds[2]-invInternal$refEnd+1)
invInternal <- invInternal[order(invInternal$newRefEnd),]
write.table(invInternal[,c("refChr","newRefEnd","newRefStart","contig","start","end","strand")],file="chr3R_inversion.tsv",sep="\t",quote=F,row.names=F,col.names=F)




inv3R <- c("NODE_99551_length_324107_cov_18.408346","NODE_79583_length_340455_cov_18.334728")
inv3R_size <- c(324107,79583)
inv3R_color <- c("red","blue")
inv3Rmapping <- mapped[mapped$contig %in% inv3R & mapped$refChr=="3R",]

pdf("chr3Rinversion.pdf",width=5,height=3)
par(mar=c(3,3,1,.5),cex.axis=.68,cex.lab=.68)
plot(0,0,col="transparent",xlim=range(inv3Rmapping$refStart,inv3Rmapping$refEnd),ylim=c(1,324107+340455),xlab="dmel",ylab="dsim",xaxs="i",yaxs="i")
abline(h=324107)

abline(v=inv3R_bounds)
arrows(inv3R_bounds[1],324107,inv3R_bounds[2],324107,lwd=5,xpd=T,lend=1)
segments(invInternal$refStart,324107,invInternal$refEnd,324107,lwd=1,col="gray68",lend=1)

padding <- 0
for (i in 1:2) {
   contigHits <- inv3Rmapping[inv3Rmapping$contig==inv3R[i],]

   for (c in 1:nrow(contigHits)) {
      rect(contigHits[c,]$refStart,1,contigHits[c,]$refEnd,324107+340455,col=rgb(t(col2rgb(inv3R_color[i]))/255,alpha=.5),bor="transparent",xpd=T)
      if (contigHits[c,]$strand=="-") {
         arrows(contigHits[c,]$refStart,contigHits[c,]$end+padding,contigHits[c,]$refEnd,contigHits[c,]$start+padding,lwd=3,col=inv3R_color[i],xpd=T)
      } else {
         arrows(contigHits[c,]$refStart,contigHits[c,]$start+padding,contigHits[c,]$refEnd,contigHits[c,]$end+padding,lwd=3,col=inv3R_color[i],xpd=T)
      }
   }

   padding <- padding + inv3R_size[i]
}
dev.off()


