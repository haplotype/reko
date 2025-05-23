# ReKo
ReKo (Reflection Knockoff) takes an input of feature matrix and produce multiple copies of knockoffs. The construction is based on a novel method, and the following paragarph is an exerpt from the manuscript in preparation. 

> In this paper, we introduce a novel method to estimate S by first constructing mirror images Y of the original feature X via Householder reflection. Their difference, W= 1/2 (X-Y), is orthogonal to X and thus T = cov(X, W) is a diagonal matrix.  We define S = \alpha T, and estimate the scalar \alpha by examining the largest eigenvalues of the scaled W.  Such an S=\hat\alpha T leads to much more powerful knockoffs, and its computation requires no optimization.  We call the knockoffs constructed via Householder reflection **reflection knockoff** (ReKo).

The basics of knockoff filter can be found here: https://web.stanford.edu/group/candes/knockoffs/.

## Usage: 
```
Usage:   reko -f filename -c number -t thread -o output_prefix [-b:c:f:o:r:s:t:]
Options: 
         -b int        number of PCs for rSVD in approximate EB-estimates [0]
         -c int        copies of knockoffs to be generated [1]
         -f str        input file contain features nxp
         -o str        output prefix [out]
         -r int        random seedL
         -s char       e: empirical Baeys estimates; p: inv_sigma = sqrt(p) [e]
            flt        inv_sigma where b~MNV(0,V) and V = inv_sigma I_p 
         -t int        number of threads [8]

Bug report: Yongtao Guan <ytguan@gmail.com>
```
