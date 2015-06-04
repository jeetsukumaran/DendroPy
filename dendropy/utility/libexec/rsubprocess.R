#!/usr/bin/env Rscript
f <- file("stdin")
open(f)
while(length(input<-readLines(f))> 0) {
    eval(parse(text=input))
}
