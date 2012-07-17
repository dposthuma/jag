
generate_qq <- function(pvalues, group, filename, adjusted)
{
	pdf.title <- paste(filename)
	pdf(file=pdf.title)
	
	plot_qq(pvalues,group)
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
#loading R file wit all functions
source(paste(r_path,"global_functions.R",sep=""))

input <- read.delim(assocfile, sep = '', header = TRUE, colClasses = "character")

if(adjusted == 'True')
{
	pvalues <- input$GC
}

if(adjusted == 'False')
{
	pvalues <- input$P
}


generate_qq(pvalues, header, filename)
