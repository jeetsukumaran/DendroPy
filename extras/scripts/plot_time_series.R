postscript("ts.ps");
ts <- read.table("t", header=TRUE)
plot(ts$iter, ts$lnL, pch=".", col="red");
hasBIter <- ts$iter[!is.na(ts$bestLnL)]
hasBLnL <- ts$bestLnL[!is.na(ts$bestLnL)]
lines(hasBIter, hasBLnL, xlab="iter", ylab="lnL", lty=1);
dev.off();
