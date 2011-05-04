### 04.19.11
# qsub -l 1hr -cwd -N prepStitching "R CMD BATCH prepStitchCoords.R"
min_match <- 500
min_contig_length <- 1000

overlapCheck <- function(focal,test) {
   test$start  <- as.numeric(as.vector(test$start))
   test$end    <- as.numeric(as.vector(test$end))
   focal$start <- as.numeric(as.vector(focal$start))
   focal$end   <- as.numeric(as.vector(focal$end))

   if ((test$start > test$end) | (focal$start > focal$end)) { FALSE }
   else if (focal$start == focal$end) { FALSE }
   else if (((test$start <= focal$start) & (test$end >= focal$end)) |
       ((test$start >= focal$start) & (test$end <= focal$end)) |
       ((test$start <= focal$start) & (test$end >= focal$start)) |
       ((test$start <= focal$end)   & (test$end >= focal$end))) { TRUE } 
   else { FALSE }
}


### READ IN BLAT DATA
####################################################################################################
blat2dmel <- read.csv("contigs_63_Dmel.psl",header=F,sep="\t",skip=6,as.is=T)
names(blat2dmel) <- c("match","mismatch","repMatch","Ns","gapCount","gapBases","ref_gapCount","ref_gapBases","strand","contig","contigLength","start","end","refChr","refSize","refStart","refEnd","blocks","blockSizes","blockStarts","ref_blockStarts")

blat2dmel$alnLength <- blat2dmel$match + blat2dmel$mismatch + blat2dmel$repMatch + blat2dmel$Ns

blat2dsim <- read.csv("contigs_63_Dsim.psl",header=F,sep="\t",skip=6,as.is=T)
names(blat2dsim) <- c("match","mismatch","repMatch","Ns","gapCount","gapBases","ref_gapCount","ref_gapBases","strand","contig","contigLength","start","end","refChr","refSize","refStart","refEnd","blocks","blockSizes","blockStarts","ref_blockStarts")
blat2dsim$alnLength <- blat2dsim$match + blat2dsim$mismatch + blat2dsim$repMatch + blat2dsim$Ns


### FILTER DATA
####################################################################################################
# [1] "Dmel starting with: 202117"
cat(sprintf("Dmel starting with: %d",nrow(blat2dmel)))

# [1] "Dmel removing <1000 bp hits: 202117"
# [1] "Dmel retaining: 28626"
bad <- blat2dmel$alnLength < min_match
blat2dmel <- blat2dmel[!bad,]
cat(sprintf("Dmel removing <%d bp hits: %d",min_match,length(bad)))
cat(sprintf("Dmel retaining: %d",nrow(blat2dmel)))
rm(bad)

### FILTER CONTIGS FOR STITCHING (contigs > 10000)
####################################################################################################
# [1] "Dmel|Dsim contigs >= 10k starting with: 1435"
# [1] "Dmel&Dsim contigs >= 10k starting with: 1426"
contigs4stitching <- unique(blat2dmel[blat2dmel$contigLength >= min_contig_length,"contig"])
cat(sprintf("contigs hitting dmel >= %d bps: %d",min_contig_length,length(contigs4stitching)))

####################################################################################################
#chroms <- c("2R")
#chroms <- c("2L","2R","3L","3R","X","4")
chroms <- c("3R","X","4")

for (chr in chroms) {
   cat(sprintf("\n*** Working on chromosome %s ***\n",chr))
   data_dmel <- blat2dmel[blat2dmel$refChr==chr & blat2dmel$contig %in% contigs4stitching,]
   data_dmel <- data_dmel[order(data_dmel$refStart),]

   pdf(sprintf("dmel_chr%s_mapping.pdf",chr),height=3,width=5)
   par(mar=c(5,10,2,1),cex.lab=.68,cex.axis=.68)

   overlapResults <- {}
   for (i in 1:nrow(data_dmel)) {
      focal_dmel <- data_dmel[i,]
      cat(sprintf("\n%s %d..%d ::: %s:%d..%d (%d)\n", focal_dmel$contig,focal_dmel$start,focal_dmel$end,focal_dmel$refChr,focal_dmel$refStart,focal_dmel$refEnd,focal_dmel$alnLength))

      status_note <- 'not in contigs4stitching';
      if (focal_dmel$contig %in% contigs4stitching) {

         ### does this region of the contig map to a BETTER location in dmel?
         better4Contig <- blat2dmel[
            (blat2dmel$contig == focal_dmel$contig & (blat2dmel$start <= focal_dmel$start) & (blat2dmel$end >= focal_dmel$end)) &
            (blat2dmel$alnLength > focal_dmel$alnLength) &
            !(blat2dmel$refChr == focal_dmel$refChr & blat2dmel$refStart == focal_dmel$refStart & blat2dmel$refEnd == focal_dmel$refEnd),]

         overlap2Contig <- blat2dmel[
            (blat2dmel$contig == focal_dmel$contig & 
            (((blat2dmel$start <= focal_dmel$start) & (blat2dmel$end >= focal_dmel$end)) |
             ((blat2dmel$start >= focal_dmel$start) & (blat2dmel$end <= focal_dmel$end)) |
             ((blat2dmel$start <= focal_dmel$start) & (blat2dmel$end >= focal_dmel$start)) |
             ((blat2dmel$start <= focal_dmel$end)   & (blat2dmel$end >= focal_dmel$end)))) &
             (blat2dmel$alnLength > focal_dmel$alnLength) &
            !(blat2dmel$refChr == focal_dmel$refChr & blat2dmel$refStart == focal_dmel$refStart & blat2dmel$refEnd == focal_dmel$refEnd),]

         ### does this region of dmel map to a BETTER contig?
         better4dmel <- blat2dmel[
            (blat2dmel$refChr == focal_dmel$refChr & 
            (blat2dmel$refStart <= focal_dmel$refStart) & (blat2dmel$refEnd >= focal_dmel$refEnd)) &
#            (((blat2dmel$refStart <= focal_dmel$refStart) & (blat2dmel$refEnd >= focal_dmel$refEnd)) |
#             ((blat2dmel$refStart >= focal_dmel$refStart) & (blat2dmel$refEnd <= focal_dmel$refEnd)) |
#             ((blat2dmel$refStart <= focal_dmel$refStart) & (blat2dmel$refEnd >= focal_dmel$refStart)) |
#             ((blat2dmel$refStart <= focal_dmel$refEnd)   & (blat2dmel$refEnd >= focal_dmel$refEnd)))) &
             (blat2dmel$alnLength > focal_dmel$alnLength) &
            !(blat2dmel$contig == focal_dmel$contig & blat2dmel$start == focal_dmel$start & blat2dmel$end == focal_dmel$end),]

         cat("Better hits for contig\n")
         cat(sprintf("\t%s:%d-%d (%d)\n",
            better4Contig$refChr,better4Contig$refStart,better4Contig$refEnd,better4Contig$alnLength))

         cat("Better hits for dmel region\n")
         cat(sprintf("\t%s:%d-%d (%d)\n",
            better4dmel$contig,better4dmel$start,better4dmel$end,better4dmel$alnLength))

         status_note <- 'lost to better pairings';
         if ((nrow(better4Contig) == 0) & (nrow(better4dmel) == 0)) {
            status_note <- 'success';
            overlapResults <- rbind(overlapResults, focal_dmel[,c("refChr","refStart","refEnd","contig","start","end","strand")])

            data_dsim <- blat2dsim[blat2dsim$contig==focal_dmel$contig,]
            if (nrow(data_dsim)>0) {
               data_dsim$overlap <- unlist(apply(data_dsim, 1, function(y) {
                  x <- as.data.frame(t(y))
                  test_start <- as.numeric(as.vector(x$start))
                  test_end   <- as.numeric(as.vector(x$end))
                  if ((test_start <= focal_dmel$start) & (test_end >= focal_dmel$end)) { focal_dmel$end - focal_dmel$start }
                  else if ((test_start >= focal_dmel$start) & (test_end <= focal_dmel$end)) { test_end - test_start }
                  else if ((test_start <= focal_dmel$start) & (test_end-1 >= focal_dmel$start)) { test_end - focal_dmel$start }
                  else if ((test_start <= focal_dmel$end-1) & (test_end >= focal_dmel$end)) { focal_dmel$end - test_start }
                  else { 0 } }))

               data_dsim <- data_dsim[data_dsim$overlap > 0,]
#               a <- as.data.frame(table(sprintf("%s.%d.%d.%d.%d",data_dsim$refChr,data_dsim$refStart,data_dsim$refEnd,data_dsim$match,data_dsim$mismatch)))
#               ### remove these duplicates, meant for plotting dsim hits anyway
#               if (nrow(a[a$Freq > 1,])>0) { 
#                  print(a[a$Freq > 1,])
#
#                  dups <- as.data.frame(matrix(unlist(lapply(as.vector(a[a$Freq > 1,]$Var1),strsplit,split="\\.")),byrow=T,ncol=5))
#                  data_dsim$dup <- FALSE
#                  data_dsim$dup <- unlist(apply(data_dsim,1,function(y) { 
#                                       x <- as.data.frame(t(y));
#                                       print(x)
#                                       any(unlist(apply(dups,1,function(w) {
#                                          v <- as.data.frame(t(w));
#                                          print(v)
#                                          if (x$refChr==v$V1 & x$refStart==v$V2 & x$refEnd==v$V3 & x$match==v$V4 & x$mismatch==v$V5) { TRUE }
#                                          else { FALSE }
#                                       })))
#                                    }))
#
#                  as.data.frame(matrix(unlist(lapply(as.vector(a[1:5,]$Var1),strsplit,split="\\.")),byrow=T,ncol=5))
#                  data_dsim <- data_dsim[data_dsim$
#                  #dev.off()
#               }
#               row.names(data_dsim) <- sprintf("%s.%d-%d.%d.%d",data_dsim$refChr,data_dsim$refStart,data_dsim$refEnd,data_dsim$match,data_dsim$mismatch)

               data_dsim$status <- 1

               hits2dmel <- nrow(overlapsINdmel)
               hits2dsim <- nrow(data_dsim)

               plot(c(1,focal_dmel$contigLength),c(1,1),type="l",lwd=1,lend=1,xlim=c(1,focal_dmel$contigLength),ylim=c(1,hits2dmel+hits2dsim+1),
                  main="",xlab=focal_dmel$contig,ylab="",axes=F)
               rect(focal_dmel$start,1,focal_dmel$end,hits2dmel+hits2dsim+1,col="gray88",bor="transparent")
               segments(focal_dmel$start,1,focal_dmel$end,1,lwd=5,lend=1)
         
               axis(1)
               mtext(side=3,line=.05,sprintf("%s:%d..%d (%s) %d bps",focal_dmel$refChr,focal_dmel$refStart,focal_dmel$refEnd,
                  focal_dmel$strand,focal_dmel$refEnd-focal_dmel$refStart),font=2)
         
               if (hits2dmel>=1) {
                  segments(overlapsINdmel$start,c(2:hits2dmel),overlapsINdmel$end,c(2:hits2dmel),col=1:(hits2dmel-1),lwd=1,lend=1) 
                  mtext(sprintf("%s:%d (%d)",overlapsINdmel$refChr,overlapsINdmel$refStart,overlapsINdmel$refEnd-overlapsINdmel$refStart),side=2,line=.05,at=c(2:hits2dmel),cex=.68,las=2)
               }
   
               segments(data_dsim$start,c(1:hits2dsim)+hits2dmel,data_dsim$end,c(1:hits2dsim)+hits2dmel,col=1:hits2dsim,lwd=1,lend=1) 
               mtext(sprintf("%s:%d (%d)",data_dsim$refChr,data_dsim$refStart,data_dsim$overlap),side=2,line=.05,at=c(1:hits2dsim)+hits2dmel,cex=.68,las=2,font=data_dsim$status)
               abline(h=cnv_dmel+.5,col="gray68",lty=2)
            }

#         } else {
#            cat("Beaten out by\n")
#            cat(sprintf("\t%s:%d-%d (%d)\n",
#               better_than_dmel$refChr,better_than_dmel$refStart,better_than_dmel$refEnd,better_than_dmel$alnLength))
         }
      }

      cat(sprintf("\t%s\n", status_note))
   }

   write.table(as.matrix(overlapResults),file=sprintf("dmel_chr%s_mapping.tsv",chr),row.names=F,col.names=F,sep="\t",quote=F)
   dev.off()
}
