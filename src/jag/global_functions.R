plot_qq<- function(pvalues,group){
	#	par(mfrow = c(5,4))
	obs <- as.numeric(pvalues)
	#remove NA (from last line etc.)
	obs<-obs[!is.na(obs)]
	sort.pval.obs <- sort(obs, decreasing=FALSE)
	log10.pval.obs <- -log10(sort.pval.obs)
      
	y = log10.pval.obs
	
	# EXPECTED #
	L <- length(obs)
	j <- c(1:L)
	x <- -log10(j/(L+1))	
		
	plot(x, y, xlim= c(0,8), ylim=c(0,8) , xlab = "Expected -log10(P)", ylab = "Observed -log10(P)", pch=20, main=group, cex.main=1, cex.lab=1, col = "black")
	abline(0,1, col = "red", lwd=1.5)
}