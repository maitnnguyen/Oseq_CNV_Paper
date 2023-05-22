#!/usr/bin/env Rscript

# 1 list of normal normalized coverage file from S2; 2 normDB name
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args) == 0) {
  stop("Check input, need at least 1 input: list of normal coverage files", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = "normDB.rds"
}

suppressWarnings(suppressMessages({
	library(PureCN)
	library(dplyr)
	library(stringr)
	}))

cov.files = read.delim(args[1])

if (is.na(args[3])){
	cov.files$File <- as.character(cov.files$File)
	normal.files <- cov.files$File[grepl('BDNA', cov.files$File)] #cov.files$File[grep('BDNA',cov.files$File)]
} else {
	ctdnaInfo = read.delim(args[3]) %>%
	mutate(type = as.character(type),
         TP53.freq = as.numeric(TP53.freq),
         SampleID = sub('_O','_TPL',sub('_DNA.*','',sample)),
         days_min.sample1.oper1. = as.numeric(days_min.sample1.oper1.))
	ctdnaInfo$SampleID[grepl('r1Asc',ctdnaInfo$sample)] = ctdnaInfo$samplename_old[grepl('r1Asc',ctdnaInfo$sample)]
	ctdnaInfo$TP53.freq[is.na(ctdnaInfo$TP53.freq)] = 0

	zeroTP53.cov.files <- cov.files %>%
		mutate(SampleID = sub('_O','_TPL',Key),
				File = as.character(File)) %>%
		left_join(ctdnaInfo[,c('SampleID', 'TP53.freq','type')]) %>%
		filter(type %in% c('plasma','tissue','blood'),
				TP53.freq == 0)
	normal.files <- na.omit(zeroTP53.cov.files$File)
}

normDB <- createNormalDatabase(normal.files)
saveRDS(normDB, file = args[2])