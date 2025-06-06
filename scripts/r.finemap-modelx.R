args <- commandArgs(trailingOnly = TRUE) 
# Parse the arguments (they are strings by default, so you may need to convert them)
gene <- as.character(args[1])  # First argument

print(gene);

library(knockoff);
fx= paste0("geno.",gene,".tr");
xx=read.table(fx); 
xx=scale(xx); 

fnph=paste0("res.",gene,".residual"); 
yy=scan(fnph); 

for(t in c(1:2)) 
{
    pref=paste0(gene,".",t); 
    for (j in c(1:10)) {
	ko=create.second_order(xx,method="sdp");
        tx=cbind(xx,ko);
        rs=paste0("rs",seq(ncol(tx)));
        aA=rep("G", ncol(tx));
        aB=rep("T", ncol(tx));
        xout = cbind(rs, aA, aB, t(round(tx,4)));
        fngt=paste0("modelx-mgt.",pref,".",j);
        write.table(xout, fngt, row.name=F,col.name=F,quote=F);

       pout=paste0("modelx.",pref,'.',j);
       rand = ceiling(runif(1) * 1e6);
       cmd = paste0("~/work/fastBVSR/fastBVSR-linux -g ",fngt," -p ",fnph," -w 100000 -smax 100 -s 1000000 -r ", rand," -o ",pout, " > /dev/null 2>&1 &" );
        print(cmd);
        system(cmd);

        w2=stat.glmnet_coefdiff(xx,ko,yy)
        mx=cbind(mx,w2);
    }
    fout=paste0("modelx-glmnet.",pref);
    write.table(mx,fout,row.names=F,col.names=F,quote=F)
}

