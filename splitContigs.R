#!/usr/bin/env Rscript
### 04.19.11
# qsub -l 1hr -cwd -N prepStitching "R CMD BATCH prepStitchCoords.R"

source("~/R/getopts.R")
contigSplit <- opts$c
refChr <- opts$r
refStart <- opts$s
refEnd <- opts$e

### READ IN BLAT DATA
####################################################################################################
blat2dmel <- read.csv("min1000_dmelARMS_3Rinv2contigs.psl",header=F,sep="\t",skip=6,as.is=T)
names(blat2dmel) <- c("match","mismatch","repMatch","Ns","ref_gapCount","ref_gapBases","gapCount","gapBases","strand",
                      "refChr","refSize","refStart","refEnd",
                      "contig","contigLength","start","end",
                      "blocks","blockSizes","ref_blockStarts","blockStarts")

contigData <- blat2dmel[blat2dmel$contig==contigSplit & blat2dmel$refChr==refChr & blat2dmel$refStart==refStart & blat2dmel$refEnd==refEnd,]
ref_blockStarts <- as.numeric(unlist(strsplit(contigData$ref_blockStarts,split=",")))
blockStarts <- as.numeric(unlist(strsplit(contigData$blockStarts,split=",")))
if (contigData$strand == "-") { blockStarts <- contigData$contigLength - blockStarts }
blockSizes <- as.numeric(unlist(strsplit(contigData$blockSizes,split=",")))
dist2last <- ref_blockStarts[-1] - ref_blockStarts[-length(ref_blockStarts)]
where2split <- which(dist2last==max(dist2last))
#ref_blockStarts[1:where2split]
#ref_blockStarts[-(1:where2split)]

cat(sprintf("%s\t%d\t%d\t%s\t%d\t%d\t%s\n",
   refChr,ref_blockStarts[1],ref_blockStarts[where2split]+blockSizes[where2split],
   contigSplit,blockStarts[1],blockStarts[where2split]+blockSizes[where2split],contigData$strand))
cat(sprintf("%s\t%d\t%d\t%s\t%d\t%d\t%s\n",
   refChr,ref_blockStarts[where2split+1],ref_blockStarts[length(ref_blockStarts)]+blockSizes[length(ref_blockStarts)],
   contigSplit,blockStarts[where2split+1],blockStarts[length(ref_blockStarts)]+blockSizes[length(ref_blockStarts) ],contigData$strand))

