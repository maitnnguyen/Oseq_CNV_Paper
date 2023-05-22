#!/usr/bin/env Rscript

# 1 coverage file of cancer sample; 2 vcf file from Mutect; 3 normDB; 
# 4 interval.file; 5 LogR output; 6 segmentation output
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args) < 4) {
  stop("Check input, need at least 4 inputs: 1-coverage.file of cancer sample; 2-vcf.file from Mutect; 3-normDB; 4-interval.file ", call.=FALSE)
} else if (length(args)==5) {
  # default output file
  args[5] = "logR.csv"
  args[6] = "segment.csv"
  #args[7] = "segment_by_gene.csv"
  args[7] = "result.rds"
}

suppressWarnings(suppressMessages({
	library(PureCN)
	library(dplyr)
	library(stringr)
  library(dplyr)
  library(copynumber)
  library(GenomicRanges)
  library(magrittr)
	}))
normDB = readRDS(args[3])
pool = calculateTangentNormal(args[1],normDB)

#skip_to_next <- FALSE

ret = runAbsoluteCN(calculateTangentNormal(args[1],normDB),
                    tumor.coverage.file=args[1],vcf.file=args[2],
                    genome='hg38', #sampleid=param3,
                    post.optimize=T,
                    interval.file=args[4],
                    normalDB=normDB,plot.cnv=F,verbose=F,
                    sex = 'F',
		                min.ploidy = 1,
                    max.ploidy = 6,
                    test.num.copy = 0:20,
                    test.purity = seq(0.01, .99, by = 0.01)
                    )

saveRDS(ret,file = args[7])

ll = c()
for (i in 1:length(ret$results)){
    ll=c(ll,ret$results[[i]]$log.likelihood)
}
sel = which.max(ll)



# output from ret
seg=cbind.data.frame(ret$results[[sel]]$seg,ret$results[[sel]]$ML.Subclonal,
                    ret$results[[sel]]$ML.C)

############### adjusted logR (Kari's discussion)
# R(x) = (aq(x)+2(1-a))/D
# D = aT + 2(1-a)
# q(x) = DR(x)/a - 2(1-a)/a
# R'(x) = q(x)/T = R(x)/a - 2(1-a)/aT
# where:
# R(x) = raw (observed) coverage ratio (PureCN's seg.mean, gene.mean, etc.)
# R'(x) = adjusted coverage ratio (in tumor cells)
# q(x) = integer copy number in cancer cells
# D = average ploidy across all cells of tumor (of sample)
# a = sample purity
# T = tumor ploidy
# drive formula for R'(x)
# R'(x) = q(x)/T = DR(x)/aT - 2(1-a)/aT = (aT + 2(1-a))R(x)/aT - 2(1-a)/aT
#       = R(x) + 2(1-a)R(x)/aT - 2(1-a)/aT
#       = [aTR(x) + 2(1-a)R(x) - 2(1-a)]/aT
rds <- ret$results[[sel]]
R <- 2^(seg$seg.mean)

seg$adj.seg = (rds$purity*rds$ploidy*R + 2*(1-rds$purity)*(R - 1))/(rds$purity*rds$ploidy)

# prediction for subclonal
mut.pred=predictSomatic(ret)

write.table(pool,file=args[5],col.names=T,row.names=F,quote=F,sep='\t')

write.table(seg,file=args[6],col.names=T,row.names=F,quote=F,sep='\t')




