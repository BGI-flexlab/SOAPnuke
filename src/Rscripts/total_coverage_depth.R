#!/usr/bin env Rscript
args <- commandArgs(TRUE)
n<-length(args)
if (n!=2) stop("\nusage: data picture.png\nAuthor:lanzhzh\nE-mail:lanzhangzhang@genomics.org.cn")
a <- read.delim(args[1],header=TRUE)
attach(a)
a <- a[order(AverageDepth),]
a[,1] <- sort(Order)
wth <- length(a[,1])/2
if ( wth < 10 ) {
    wth <- 10
}
bitmap(file=args[2],width = wth,height =10,type = "png16")
par(ann=F,mar=c(4,4,4,4),xpd=TRUE)
plot(x=a[,1],y=a[,3],type="o",col="blue",xaxt="n",ylim=c(0,1.1*max(a[,3])),cex.axis=0.8,las=1,pch=20,xaxs="i",yaxs="i")
axis(1,at=a[,1],lab=FALSE)
text(x=a[,1],y=-5,a[,2],font=2,srt=45,cex=0.9)
title(main="Sample Quality Statistic",font.main=4)
title(xlab="Sample",ylab="Depth")
m<-max(a[,3])
n<-(m/100)*1.04
lines(x=a[,1],y=a[,4]*n,col="green")
points(x=a[,1],y=a[,4]*n,pch=15,col="green",cex=0.5)
lines(x=a[,1],y=a[,5]*n,col="red")
points(x=a[,1],y=a[,5]*n,pch=18,col="red",cex=0.8)
abline(v=a[,1],lty =2,col="gray50",xpd=FALSE)
axis(4,at=seq(0,100*n,length=6),lab=seq(0,100,by=20),cex.axis=0.8,las=1)
mtext("rate(%)",4,3)
legend("topleft",pch=c(20,15,18),col=c("blue","green","red"),legend=c("Mean depth","Depth>=20(%)","Coverage(%)"),horiz="true",cex=0.8,bty="n",lty=1,text.width=2.3)
box()

