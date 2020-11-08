################ Bernstein polynomial ###################################
#' @export
plot_b <- function(d1, d2)
{
unm1 <-  umd(sort(d1),0,max(d1)+sd(d1)/length(d1))
xvec=seq(0,10,by=0.001)
l1 <- unm1$dumd(xvec)
plot(xvec,l1, xlab=" ", col='red', ylab=" ",pch='.',xaxt='n', yaxt='n', cex.lab=2, lwd=2, cex.main=1.5, cex.sub=2)
par(new=TRUE)
unm2 <-  umd(sort(d2),0,max(d2)+sd(d2)/length(d2))
l2 <- unm2$dumd(xvec)
plot(xvec, l2,  pch='.', xlab="log(net MFI) ", col='blue', ylab=" ",  cex.lab=2, lwd=2,cex.axis=1.5, cex.main=1.5, cex.sub=2)
legend("bottomright",c("HVTN 097","HVTN 100"),lwd = rep(2,2),lty=rep(1,2),
       col=c('red','blue'),pt.cex=1.4, cex=1.4, bty='n')
title(main="Turnbull and Ghosh (2014)'s  density estimators")
grid(lty=1)
}
