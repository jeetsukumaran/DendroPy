postscript("ts.ps");
ts <- read.table("t", header=TRUE)
plot(ts$iter, ts$lnL, pch=".", col="red", xlab="# trees scored", ylab="lnL", cex.lab=1.5);
hasBIter <- ts$iter[!is.na(ts$bestLnL)]
hasBLnL <- ts$bestLnL[!is.na(ts$bestLnL)]
lines(hasBIter, hasBLnL, lty=1);
dev.off();
