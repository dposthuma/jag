dist_plot<-function(emp_file,outfile,prefix,pheno){
  #################################### DISTRIBUTION PLOTS  #####################################
  ## Met de onderstaande code kan je voor iedere groep die je hebt getest de distributie van de permutaties plotten.
  ## Ik heb deze code aangepast op de JAG output files
  ##  Ik heb de code getest op 100 permutaties. Input data staat in 'input_output_plots'. 

  ## READ IN FILE WITH EMP PVALUES AND GET ORDER FROM LOWEST Pemp to HIGHEST Pemp
  distr.sumlogs <- read.table(emp_file, sep="\t", header=TRUE, colClasses='character')
  ordered.sumlogs <- order(as.numeric(distr.sumlogs[,4]))

  ## GET PERMUTED SUMLOGS (NOTE: de naamgeving van deze file is nu variabel!)
  permuted.sumlogs <- read.delim(outfile, sep = "\t", colClasses ='character') 

  nr.groups <- ncol(permuted.sumlogs)-1
  pdf.title <- paste(prefix,"distribution_sumlogs.P",pheno,".pdf", sep = '')
  pdf(file=pdf.title)
  
  for(i in 1:length(ordered.sumlogs))
  {	
	  # OPEN PNG > note that name of file is not specific now!

	  print(paste("plotting group ", distr.sumlogs[ordered.sumlogs[i],1], sep=''))	
	  
	  # Generate 100 (or another number) equal bins and take the min and max of permuted data as extremes of bins  
	  seq.data <- seq(from = min(as.numeric(permuted.sumlogs[,ordered.sumlogs[i]])),to = max(as.numeric(permuted.sumlogs[,ordered.sumlogs[i]])), length.out = 100)
	  
	  values <- c(as.numeric(permuted.sumlogs[,ordered.sumlogs[i]]),distr.sumlogs[ordered.sumlogs[i],2])
	  min <- min(as.numeric(values))/100*95
	  max <- max(as.numeric(values))/100*105
	  
	  # plot histogram of permuted distribution
	  hist(as.numeric(permuted.sumlogs[,ordered.sumlogs[i]]), seq.data, main = distr.sumlogs[ordered.sumlogs[i],1], xlab = "sumlog values", xlim = c(min, max), cex.main=1.8, cex.lab=1.5, col="lightgrey")
	  abline(v =distr.sumlogs[ordered.sumlogs[i],2], col="red", lty=2, cex=2, lwd=2)  
  }
  dev.off()
}

#read all arguments given by python
#DO NOT REMOVE!
args=(commandArgs(TRUE))
if(length(args)==0){
    print("No arguments supplied.")
}else{
    for(i in 1:length(args)){
         eval(parse(text=args[[i]]))
    }
}
dist_plot(emp_file,out_file,prefix,pheno)