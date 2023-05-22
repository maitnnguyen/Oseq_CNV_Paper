#!/usr/bin/env Rscript

# 1 coverage file from S1, 2 interval.file, 3 output
# interval.file here has to contain GC and mappability value
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args) < 2) {
  stop("Check input, need at least 2 inputs: raw coverage file, interval.file, and output (optional)", call.=FALSE)
} else if (length(args)==2) {
  # default output file
  args[3] = "out.txt"
}

suppressWarnings(suppressMessages({
	library(PureCN)
	library(dplyr)
	library(stringr)
}))
# correct GC content and mappability
correctCoverageBias(coverage.file=args[1] , interval.file = args[2], 
                      output.file = args[3])
  