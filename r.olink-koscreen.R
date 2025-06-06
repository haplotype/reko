xx=read.table("olink.14752.residual"); 
xx=as.matrix(xx); 
##xx is n by p; 
fnph="age.14752.residual"; 
yy=scan(fnph);
##yy is n-vector;
library(KnockoffScreen); 
library(knockoff); 
    
mk=create.MK(xx,seq(ncol(xx)),10); 
for(j in c(1:10)) {
    ko=as.matrix(mk[[j]]); 

    tx=cbind(xx,ko);
    rs=paste0("rs",seq(ncol(tx)));
    aA=rep("G", ncol(tx));
    aB=rep("T", ncol(tx));
    xout = cbind(rs, aA, aB, t(round(tx,4)));
    fngt=paste0("koscreen-mgt.",j,'.txt');
    write.table(xout, fngt, row.name=F,col.name=F,quote=F);

   pout=paste0("koscreen.",j); 
   rand = ceiling(runif(1) * 1e6); 
   cmd = paste0("~/work/fastBVSR/fastBVSR-linux -g ",fngt," -p ",fnph," -w 100000 -s 1000000 -r ", rand," -o ",pout, " > /dev/null 2>&1 &" ); 
    print(cmd); 
    system(cmd); 

    w2=stat.glmnet_coefdiff(xx,ko,yy); 
    fn=paste0("koscreen.glmnet.",j); 
    write.table(w2, fn, row.names=F,col.names=F,quote=F); 
}
