library(knockoff); 
#cmd="~/work/knockoff/reko -f olink.14752.residual -c 10 -o olink.14572";
#system(cmd); 
##this is to call a C program reko to construct 10 copies of reko knockoffs. 


xx=read.table("olink.14752.residual"); 
##xx is an  n by p feature matrix; 
ph="age.14752.residual"; 
yy=scan(ph);
##yy is an n-vector phenotype;
res=c(); 
for(t in c(1:10)) {
    kk=read.table(paste0("olink.14752.res.",t,".ko"));  
    ## read reko knockoffs from files.  
    ## you may construct knockoff with existing method using 
    ##  kk = create.second_order(xx,method="sdp"); 

    tx=cbind(xx,kk);
    rs=paste0("rs",seq(ncol(tx)));
    aA=rep("G", ncol(tx));
    aB=rep("T", ncol(tx));
    xout = cbind(rs, aA, aB, t(round(tx,4)));
    gt=paste0("olink.mgt-ko.",t,".txt");
    write.table(xout, gt, row.name=F,col.name=F,quote=F);
    ##prepare input file for fastBVSR; 

    pout=paste0("olink-bvsr.",t);
    rand=ceiling(runif(1)*1e6); 
    cmd = paste0("~/work/fastBVSR/fastBVSR-linux -g ",gt," -p ",ph," -r ",rand, " -w 100000 -s 1000000 --pi-max 0.9 -o ",pout, " > /dev/null 2>&1 &" );
    print(cmd);
    system(cmd);
    ##run fastBVSR, the output of interest is *.beta.txt (see below);
    ##fastBVSR is available at https://github.com/zhouquan34/fastBVSR 

    w2=stat.glmnet_coefdiff(xx,kk,yy);
    fnsu = paste0("olink-glmnet-ko.",t); 
    write(round(w2,4), fnsu, 1); 
    ##run Elastic-Net
}

stop("wait until fastBVSR finish");
library(knockoff); 
pipx=c(); 
cset=c(1:10); 
for(t in cset) {
    fn=paste0("olink-bvsr.",t,".beta.txt");
    print(fn); 
    xx=read.table(fn,1);
    p=nrow(xx)/2;    
    tp = abs(xx[1:p,4])-abs(xx[1:p+p,4]);
    pipx=cbind(pipx,tp); 
}
##pipx is a p by 10 matrix of knockoff statistics from knockoff filter ReKo:fastBVSR:PIP

w2=c(); 
for(t in cset) {
    fn=paste0("olink-glmnet-ko.",t);
    print(fn); 
    yy=scan(fn);
    w2=cbind(w2,yy); 
}
##w2 is a p by 10 matrix of knockoff statistics from knockoff filter ReKo:glmnet:Beta

mx=apply(pipx,1,mean); 
my=apply(w2,1,mean); 
for(alpha in c(0.05,0.10)) {
    thx=knockoff.threshold(as.vector(pipx),alpha);  
    thy=knockoff.threshold(as.vector(w2),alpha);  

    w=which(mx>=thx & my>=thy);
    print(c(alpha,length(w))); 
}
## w is the index of significant features from Combined-Filter at nominal FDR alpha. 
