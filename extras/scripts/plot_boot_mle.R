postscript("boot_mle.ps");
ts <- read.table("booti_ml_reps.txt", header=TRUE)
gts <- read.table("gar_ml_reps.txt", header=TRUE)

minLnL = min(gts$zer)
maxLnL = max(ts$best)
minIter = 0
maxIter = max(ts$iter, gts$iter)
plot(gts$iter, gts$zer, pch=".", col="red", xlab="# trees scored", ylab="lnL", cex.lab=1.5, ylim=c(minLnL, maxLnL), xlim=c(minIter,maxIter));
points(ts$iter, ts$zer, pch=".", col="green", xlab="# trees scored", ylab="lnL", cex.lab=1.5);

hasoneIter <- ts$iter[!is.na(ts$one)]
hasone <- ts$one[!is.na(ts$one)]
points(hasoneIter, hasone, pch=".", col="green")


hastwoIter <- ts$iter[!is.na(ts$two)]
hastwo <- ts$two[!is.na(ts$two)]
points(hastwoIter, hastwo, pch=".", col="green")

hasthrIter <- ts$iter[!is.na(ts$thr)]
hasthr <- ts$thr[!is.na(ts$thr)]
points(hasthrIter, hasthr, pch=".", col="green")

hasfouIter <- ts$iter[!is.na(ts$fou)]
hasfou <- ts$fou[!is.na(ts$fou)]
points(hasfouIter, hasfou, pch=".", col="green")

hasfivIter <- ts$iter[!is.na(ts$fiv)]
hasfiv <- ts$fiv[!is.na(ts$fiv)]
points(hasfivIter, hasfiv, pch=".", col="green")

hassixIter <- ts$iter[!is.na(ts$six)]
hassix <- ts$six[!is.na(ts$six)]
points(hassixIter, hassix, pch=".", col="green")

hassevIter <- ts$iter[!is.na(ts$sev)]
hassev <- ts$sev[!is.na(ts$sev)]
points(hassevIter, hassev, pch=".", col="green")

haseigIter <- ts$iter[!is.na(ts$eig)]
haseig <- ts$eig[!is.na(ts$eig)]
points(haseigIter, haseig, pch=".", col="green")

hasninIter <- ts$iter[!is.na(ts$nin)]
hasnin <- ts$nin[!is.na(ts$nin)]
points(hasninIter, hasnin, pch=".", col="green")


lines(ts$iter, ts$best, lty=1, lwd=2,col="blue");




hasoneIter <- gts$iter[!is.na(gts$one)]
hasone <- gts$one[!is.na(gts$one)]
points(hasoneIter, hasone, pch=".", col="red")


hastwoIter <- gts$iter[!is.na(gts$two)]
hastwo <- gts$two[!is.na(gts$two)]
points(hastwoIter, hastwo, pch=".", col="red")

hasthrIter <- gts$iter[!is.na(gts$thr)]
hasthr <- gts$thr[!is.na(gts$thr)]
points(hasthrIter, hasthr, pch=".", col="red")

hasfouIter <- gts$iter[!is.na(gts$fou)]
hasfou <- gts$fou[!is.na(gts$fou)]
points(hasfouIter, hasfou, pch=".", col="red")




hasfivIter <- gts$iter[!is.na(gts$fiv)]
hasfiv <- gts$fiv[!is.na(gts$fiv)]
points(hasfivIter, hasfiv, pch=".", col="red")

hassixIter <- gts$iter[!is.na(gts$six)]
hassix <- gts$six[!is.na(gts$six)]
points(hassixIter, hassix, pch=".", col="red")



hassevIter <- gts$iter[!is.na(gts$sev)]
hassev <- gts$sev[!is.na(gts$sev)]
points(hassevIter, hassev, pch=".", col="red")

haseigIter <- gts$iter[!is.na(gts$eig)]
haseig <- gts$eig[!is.na(gts$eig)]
points(haseigIter, haseig, pch=".", col="red")

hasninIter <- gts$iter[!is.na(gts$nin)]
hasnin <- gts$nin[!is.na(gts$nin)]
points(hasninIter, hasnin, pch=".", col="red")


lines(gts$iter, gts$best, lty=1, lwd=2,col="black");


dev.off();
