
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

raw_p_table<-read.table(p_values_file)
unique.sets <- unique(raw_p_table[,1])
pdf(file=image_filename)

for (set in 1:length(unique.sets))
{
	p_per_set=raw_p_table[raw_p_table[,1]==unique.sets[set],2]
	plot_qq(p_per_set,unique.sets[set])
}
dev.off()