#!/usr/bin/env Rscript
f <- file("stdin")
open(f)
while(length(line<-readLines(f))> 0) {
    eval(parse(text=line))
}
