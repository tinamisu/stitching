#!/usr/bin/env Rscript
### 04.19.11
# qsub -l 1hr -cwd -N prepStitching "R CMD BATCH prepStitchCoords.R"
min_match <- 500
min_contig_length <- 1000
chimera_gap <- 10000

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

plotHits <- function(focal_dmel,overlap2Contig,data_dsim,status) {
   hits2dmel <- nrow(overlap2Contig)
   hits2dsim <- nrow(data_dsim)

   if (hits2dmel > 0 | hits2dsim > 0) {
      plot(c(1,focal_dmel$contigLength),c(1,1),type="l",lwd=1,lend=1,xlim=c(1,focal_dmel$contigLength),ylim=c(1,hits2dmel+hits2dsim+1),
         main="",xlab=focal_dmel$contig,ylab="",axes=F)
      rect(focal_dmel$start,1,focal_dmel$end,hits2dmel+hits2dsim+1,col="gray88",bor="transparent")
      segments(focal_dmel$start,1,focal_dmel$end,1,lwd=5,lend=1)
      
      axis(1)
      mtext(side=3,line=.05,sprintf("%s:%d..%d (%s) %d bps\n%s",focal_dmel$refChr,focal_dmel$refStart,focal_dmel$refEnd,
         focal_dmel$strand,focal_dmel$refEnd-focal_dmel$refStart,status),font=2)
      abline(h=1+hits2dmel+.5,col="gray68",lty=2)
      
      if (hits2dmel>0) {
         segments(overlap2Contig$start,c(1:hits2dmel)+1,overlap2Contig$end,c(1:hits2dmel)+1,col=1:hits2dmel,lwd=1,lend=1) 
         mtext(sprintf("%s:%d (%d)",overlap2Contig$refChr,overlap2Contig$refStart,overlap2Contig$refEnd-overlap2Contig$refStart),side=2,line=.05,at=c(1:hits2dmel)+1,cex=.68,las=2)
      }
      
      if (hits2dsim>0) {
         segments(data_dsim$start,c(1:hits2dsim)+hits2dmel+1,data_dsim$end,c(1:hits2dsim)+hits2dmel+1,col=1:hits2dsim,lwd=1,lend=1) 
         mtext(sprintf("%s:%d (%d)",data_dsim$refChr,data_dsim$refStart,data_dsim$overlap),side=4,line=.05,at=c(1:hits2dsim)+hits2dmel,cex=.68,las=2,font=data_dsim$status)
      }
   }
}


### READ IN BLAT DATA
####################################################################################################
blat2dmel <- read.csv("min1000_dmelARMS_3Rinv2contigs.psl",header=F,sep="\t",skip=6,as.is=T)
names(blat2dmel) <- c("match","mismatch","repMatch","Ns","ref_gapCount","ref_gapBases","gapCount","gapBases","strand",
                      "refChr","refSize","refStart","refEnd",
                      "contig","contigLength","start","end",
                      "blocks","blockSizes","ref_blockStarts","blockStarts")
#blat2dmel$alnLength <- blat2dmel$match + blat2dmel$mismatch + blat2dmel$repMatch + blat2dmel$Ns

blat2dsim <- blat2dmel[1,]
#blat2dsim <- read.csv("min1000_dsimARMS2contigs.psl",header=F,sep="\t",skip=6,as.is=T)
#names(blat2dsim) <- c("match","mismatch","repMatch","Ns","ref_gapCount","ref_gapBases","gapCount","gapBases","strand",
#                      "refChr","refSize","refStart","refEnd",
#                      "contig","contigLength","start","end",
#                      "blocks","blockSizes","ref_blockStarts","blockStarts")
#blat2dsim$alnLength <- blat2dsim$match + blat2dsim$mismatch + blat2dsim$repMatch + blat2dsim$Ns

### break apart blat hits if they appear to be chimeras
blat2dmel$alnLength <- 0
blat2dmel$ref_maxGapSize <- 0
chimerasSplit <- {}
chimeras <- {}
for (i in 1:nrow(blat2dmel)) {
   blockSizes <- as.numeric(unlist(strsplit(blat2dmel[i,]$blockSizes,split=",")))
   blat2dmel[i,]$alnLength <- sum(blockSizes)

   if (blat2dmel[i,]$blocks > 1) {
      ### DETERMINE distances between blocks
      contig_blockStarts <- as.numeric(unlist(strsplit(blat2dmel[i,]$blockStarts,split=",")))
      ref_blockStarts <- as.numeric(unlist(strsplit(blat2dmel[i,]$ref_blockStarts,split=",")))

      if (blat2dmel[i,]$strand == "-") { 
         ref_blockStarts <- blat2dmel[i,]$refSize - ref_blockStarts 
         blocks <- as.data.frame(cbind(ref_blockStarts-blockSizes,ref_blockStarts,blockSizes,rev(contig_blockStarts),rev(contig_blockStarts+blockSizes)))
      } else {
         blocks <- as.data.frame(cbind(ref_blockStarts,ref_blockStarts+blockSizes,blockSizes,contig_blockStarts,contig_blockStarts+blockSizes))
      }
      names(blocks) <- c("start","end","size","contig_start","contig_end")

      blocks <- blocks[order(blocks$start),]
      blocks$dist2last <- c(0,blocks[-1,]$end - blocks[-nrow(blocks),]$start)
      blat2dmel[i,]$ref_maxGapSize <- max(blocks$dist2last) 


      ### DETERMINE split positions
      ### update start,end and refStart,refEnd positions
      where2split <- which(blocks$dist2last >= chimera_gap)
      last_w <- 1
      for (w in where2split) {
         newHit <- blat2dmel[i,]
         newHit$refStart <- blocks[last_w,]$start
         newHit$refEnd   <- blocks[w-1,]$end
         newHit$start    <- blocks[last_w,]$contig_start
         newHit$end      <- blocks[w-1,]$contig_end

         newHit$alnLength <- sum(blocks[last_w:(w-1),]$size)
         chimerasSplit <- rbind(chimerasSplit,newHit)
         last_w <- w
      }

      ### RHS of block
      if (length(where2split)>0) {
         newHit <- blat2dmel[i,]
         newHit$refStart <- blocks[last_w,]$start
         newHit$refEnd   <- blocks[nrow(blocks),]$end
         newHit$start    <- blocks[last_w,]$contig_start
         newHit$end      <- blocks[nrow(blocks),]$contig_end

         newHit$alnLength <- sum(blocks[last_w:nrow(blocks),]$size)
         chimerasSplit <- rbind(chimerasSplit,newHit)

         chimeras <- c(chimeras,i)
         cat(sprintf("(%d) %s ::: %d..%d ~ %s %d..%d\n",length(where2split),
                     blat2dmel[i,]$contig,blat2dmel[i,]$start,blat2dmel[i,]$end,
                     blat2dmel[i,]$refChr,blat2dmel[i,]$refStart,blat2dmel[i,]$refEnd))
         #print(where2split)
         #print(blocks[blocks$dist2last >= chimera_gap,])
         #cat("\n")
      }
   }
}

#chimera_gap <- 9900
#blat2dmel[blat2dmel$ref_maxGapSize>=9900 & blat2dmel$ref_maxGapSize<=10000,][1,]

pdf("chimera_cds.pdf")
plot(ecdf(blat2dmel$ref_maxGapSize))
abline(v=chimera_gap,col="red")
mtext(sprintf("#chimeras %d with gap size >= %d",nrow(blat2dmel[blat2dmel$ref_maxGapSize>=chimera_gap,]),chimera_gap),line=.05,side=3,col="blue",font=2)
qqnorm(blat2dmel$ref_maxGapSize)
abline(h=chimera_gap,col="red")
dev.off()

### remove chimeras from blat2dmel and tack on the newly split chimeras
cat(sprintf("Detected %d chimeras with a gap in blocks >= %d: broken into %d\n",length(chimeras),chimera_gap,nrow(chimerasSplit)))
blat2dmel <- rbind(blat2dmel[-chimeras,],chimerasSplit)


### FILTER DATA
####################################################################################################
# [1] "Dmel starting with: 202117"
cat(sprintf("Dmel starting with: %d\n",nrow(blat2dmel)))
cat(sprintf("Dsim starting with: %d\n",nrow(blat2dsim)))

bad <- blat2dmel$alnLength < min_match
blat2dmel <- blat2dmel[!bad,]
cat(sprintf("Dmel removing <%d bp hits: %d\n",min_match,length(bad)))
cat(sprintf("Dmel retaining: %d\n",nrow(blat2dmel)))
rm(bad)

# ### remove exact duplicates wrt contig coordinates
# tmp <- as.data.frame(blat2dmel[,c("contig","start","end")])
# dups <- which(duplicated(tmp) | duplicated(tmp,fromLast = TRUE))
# bad <- blat2dmel$alnLength < min_match
# rm(bad)

bad <- blat2dsim$alnLength < min_match
blat2dsim <- blat2dsim[!bad,]
cat(sprintf("Dsim removing <%d bp hits: %d\n",min_match,length(bad)))
cat(sprintf("Dsim retaining: %d\n",nrow(blat2dsim)))
rm(bad)


### FILTER CONTIGS FOR STITCHING (contigs > 10000)
####################################################################################################
contigs4stitching <- unique(blat2dmel[blat2dmel$contigLength >= min_contig_length,"contig"])
cat(sprintf("contigs hitting dmel >= %d bps: %d\n",min_contig_length,length(contigs4stitching)))

####################################################################################################
#chroms <- c("4")
chroms <- c("2L","2R","3L","3R","3Rinverted","X","4")
for (chr in chroms) {
   data_dmel <- blat2dmel[blat2dmel$refChr==chr & blat2dmel$contig %in% contigs4stitching,]
   data_dmel <- data_dmel[order(data_dmel$refStart),]
   cat(sprintf("\n*** Working on chromosome %s ***\nTotal hits on dmel %d\n",chr,nrow(data_dmel)))

   if (nrow(data_dmel)>0) {
      pdf(sprintf("chr%s_mapping-%d-%d.pdf",chr,min_contig_length,min_match),height=3,width=8)
      par(mar=c(5,10,2,10),cex.lab=.68,cex.axis=.68)
   
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
                ((blat2dmel$start <= focal_dmel$start) & (blat2dmel$end+1 >= focal_dmel$start)) |
                ((blat2dmel$start <= focal_dmel$end+1) & (blat2dmel$end >= focal_dmel$end)))) &
                (blat2dmel$alnLength > focal_dmel$alnLength) &
               !(blat2dmel$refChr == focal_dmel$refChr & blat2dmel$refStart == focal_dmel$refStart & blat2dmel$refEnd == focal_dmel$refEnd),]
   
            ### does this region of dmel map to a BETTER contig?
            better4dmel <- blat2dmel[
               (blat2dmel$refChr == focal_dmel$refChr & 
               (blat2dmel$refStart <= focal_dmel$refStart) & (blat2dmel$refEnd >= focal_dmel$refEnd)) &
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
                  data_dsim$status <- 1
               } 

               plotHits(focal_dmel,overlap2Contig,data_dsim,status_note)
            }
         }
         cat(sprintf("\t%s\n", status_note))
      }
   
      write.table(as.matrix(overlapResults),file=sprintf("chr%s_mapping-%d-%d.tsv",chr,min_contig_length,min_match),row.names=F,col.names=F,sep="\t",quote=F)
      dev.off()
   }
}
