
plot.overlap.l <- function(a65,a70,a75,a80) {
	a65[,3] <- log(a65[,3],base=2); a70[,3] <- log(a70[,3],base=2); a75[,3] <- log(a75[,3],base=2); a80[,3] <- log(a80[,3],base=2)
	#a65[,3] <- log(a65[,3],base=10); a70[,3] <- log(a70[,3],base=10); a75[,3] <- log(a75[,3],base=10); a80[,3] <- log(a80[,3],base=10)
	maxL <- max(a65[,2],a70[,2],a75[,2],a80[,2])
	maxC <- max(a65[,3],a70[,3],a75[,3],a80[,3])
	
	plot(1,1,col="white",xlim=c(0,maxL),ylim=c(0,maxC),xlab="Contig length",ylab="Count")
	for (i in 1:dim(a65)[1]) {
		lines(c(a65[i,1],a65[i,2]),c(a65[i,3],a65[i,3]),col="red")
	}
	for (i in 1:dim(a70)[1]) {
		lines(c(a70[i,1],a70[i,2]),c(a70[i,3],a70[i,3]),col="green")
	}
	for (i in 1:dim(a75)[1]) {
		lines(c(a75[i,1],a75[i,2]),c(a75[i,3],a75[i,3]),col="blue")
	}
	for (i in 1:dim(a80)[1]) {
		lines(c(a80[i,1],a80[i,2]),c(a80[i,3],a80[i,3]),col="orange")
	}
	
	legend("topright", c("65","70","75","80"),col=c("red","green","blue","orange"), lty=1)

}


plot.overlap.l.connected <- function(a65,a70,a75,a80) {
	a65[,3] <- log(a65[,3],base=10); a70[,3] <- log(a70[,3],base=10); a75[,3] <- log(a75[,3],base=10); a80[,3] <- log(a80[,3],base=10)
	maxL <- max(a65[,2],a70[,2],a75[,2],a80[,2])
	maxC <- max(a65[,3],a70[,3],a75[,3],a80[,3])
	
	a65x <- numeric(0); a65y <- numeric(0);
	a70x <- numeric(0); a70y <- numeric(0);
	a75x <- numeric(0); a75y <- numeric(0);
	a80x <- numeric(0); a80y <- numeric(0);
	
	plot(1,1,col="white",xlim=c(0,maxL),ylim=c(0,maxC),xlab="Contig length",ylab="Count")
	for (i in 1:dim(a65)[1]) {
		a65x <- c(a65x,a65[i,1],a65[i,2])
		a65y <- c(a65y,a65[i,3],a65[i,3])
		#lines(c(a65[i,1],a65[i,2]),c(a65[i,3],a65[i,3]),col="red")
	}
	for (i in 1:dim(a70)[1]) {
		a70x <- c(a70x,a70[i,1],a70[i,2])
		a70y <- c(a70y,a70[i,3],a70[i,3])
		#lines(c(a65[i,1],a65[i,2]),c(a65[i,3],a65[i,3]),col="red")
	}
	for (i in 1:dim(a75)[1]) {
		a75x <- c(a75x,a75[i,1],a75[i,2])
		a75y <- c(a75y,a75[i,3],a75[i,3])
		#lines(c(a65[i,1],a65[i,2]),c(a65[i,3],a65[i,3]),col="red")
	}
	for (i in 1:dim(a80)[1]) {
		a80x <- c(a80x,a80[i,1],a80[i,2])
		a80y <- c(a80y,a80[i,3],a80[i,3])
		#lines(c(a65[i,1],a65[i,2]),c(a65[i,3],a65[i,3]),col="red")
	}

	lines(a65x,a65y,col="red")
	lines(a70x,a70y,col="green")
	lines(a75x,a75y,col="blue")
	lines(a80x,a80y,col="orange")
}