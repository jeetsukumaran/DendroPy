postscript("boot_mle.ps");
ts <- read.table("ml_forr.txt", header=TRUE)
plot(ts$iter, ts$zer, pch=".", col="red", xlab="# trees scored", ylab="lnL", cex.lab=1.5);

hasoneIter <- ts$iter[!is.na(ts$one)]
hasone <- ts$one[!is.na(ts$one)]
points(hasoneIter, hasone, pch=".", col="red")


hastwoIter <- ts$iter[!is.na(ts$two)]
hastwo <- ts$two[!is.na(ts$two)]
points(hastwoIter, hastwo, pch=".", col="red")

hasthrIter <- ts$iter[!is.na(ts$thr)]
hasthr <- ts$thr[!is.na(ts$thr)]
points(hasthrIter, hasthr, pch=".", col="red")

hasfouIter <- ts$iter[!is.na(ts$fou)]
hasfou <- ts$fou[!is.na(ts$fou)]
points(hasfouIter, hasfou, pch=".", col="red")

hasfivIter <- ts$iter[!is.na(ts$fiv)]
hasfiv <- ts$fiv[!is.na(ts$fiv)]
points(hasfivIter, hasfiv, pch=".", col="red")

hassixIter <- ts$iter[!is.na(ts$six)]
hassix <- ts$six[!is.na(ts$six)]
points(hassixIter, hassix, pch=".", col="red")

hassevIter <- ts$iter[!is.na(ts$sev)]
hassev <- ts$sev[!is.na(ts$sev)]
points(hassevIter, hassev, pch=".", col="red")

haseigIter <- ts$iter[!is.na(ts$eig)]
haseig <- ts$eig[!is.na(ts$eig)]
points(haseigIter, haseig, pch=".", col="red")

hasninIter <- ts$iter[!is.na(ts$nin)]
hasnin <- ts$nin[!is.na(ts$nin)]
points(hasninIter, hasnin, pch=".", col="red")


lines(ts$iter, ts$best, lty=1, lwd=2);


dev.off();
