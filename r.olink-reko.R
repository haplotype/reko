xx=read.table("olink.14752.residual"); 
##xx is n by p; 
fnph="age.14752.residual"; 
yy=scan(fnph);
##yy is n-vector;
library(knockoff); 
    
rand1 = ceiling(runif(1) * 1e6); 
cmd = paste0("~/work/knockoff/reko2 -f olink.14752.residual -b 400 -o test -r ", rand1); 
print(cmd); 
system(cmd); 
for (j in c(1:10)) {
    fnko=paste0("test.",j,".ko"); 
    ko=read.table(fnko);

    tx=cbind(xx,ko);
    rs=paste0("rs",seq(ncol(tx)));
    aA=rep("G", ncol(tx));
    aB=rep("T", ncol(tx));
    xout = cbind(rs, aA, aB, t(round(tx,4)));
    fngt=paste0("reko-mgt.",j,'.txt');
    write.table(xout, fngt, row.name=F,col.name=F,quote=F);

   pout=paste0("reko.",j); 
   rand = ceiling(runif(1) * 1e6); 
   cmd = paste0("~/work/fastBVSR/fastBVSR-linux -g ",fngt," -p ",fnph," -w 100000 -s 1000000 -r ", rand," -o ",pout, " > /dev/null 2>&1 &" ); 
    print(cmd); 
    system(cmd); 

    w2=stat.glmnet_coefdiff(xx,ko,yy)
    fout=paste0("reko.glmnet.",j);
    write.table(w2,fout,row.names=F,col.names=F,quote=F)
}

