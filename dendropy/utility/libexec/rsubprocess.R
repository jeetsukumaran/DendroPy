#!/usr/bin/env Rscript

f <- file("stdin")
open(f)
.rsubprocess.verbose = F
while(length(line <- readLines(f,n=1)) > 0) {
    if (.rsubprocess.verbose) {
        write(paste(">", line), stderr())
    }
    eval(parse(text=line))
}
