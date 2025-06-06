args <- commandArgs(trailingOnly = TRUE) 
# Parse the arguments (they are strings by default, so you may need to convert them)
gene <- as.character(args[1])  # First argument

print(gene);

library(susieR);
library(knockoff);
fx= paste0("geno.",gene,".tr");
xx=read.table(fx); 
xx=scale(xx); 

fnph=paste0("res.",gene,".residual"); 
yy=scan(fnph); 

su=susie(xx,yy, max_iter=1000, L=40);
w2=su$pip;
fout=paste0("susie-pip.",gene);
write.table(w2,fout,row.names=F,col.names=F,quote=F)

