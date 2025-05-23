# ReKo
ReKo (Reflection Knockoff) takes an input of feature matrix and produce multiple copies of knockoffs. The construction is based on a novel method, and the following paragarph is an exerpt from the manuscript in preparation. 

> In this paper, we introduce a novel method to estimate S by first constructing mirror images Y of the original feature X via Householder reflection. Their difference, W= 1/2 (X-Y), is orthogonal to X and thus T = cov(X, W) is a diagonal matrix.  We define S = \alpha T, and estimate the scalar \alpha by examining the largest eigenvalues of the scaled W.  Such an S=\hat\alpha T leads to much more powerful knockoffs, and its computation requires no optimization.  We call the knockoffs constructed via Householder reflection **reflection knockoff** (ReKo).

The basics of knockoff filter can be found here: https://web.stanford.edu/group/candes/knockoffs/.

## Usage: 
./reko -f filename -c number -t thread -o output_prefix
1. -f str  input file contain feature matrix of n by p. 
2. -b num  number of eigenpairs for rSVD in approximate EB estimates [0]. 
3. -c num  copies of reflection knockoff to be constructed [10].
4. -t num  number of threads [8].
5. -o str  output prefix [out].
6. -r int  random seed.
7. -s chr  e: EB estimates; p: inv_sigma = sqrt(p) [e]
      flt  inv_sigma to specify b~MNV(0,V) with V = inv_sigma I_p.  
