# ReKo
ReKo (Reflection Knockoff) takes an input of feature matrix and produce multiple copies of knockoffs. The construction is based on a novel method, and the following paragarph is an exerpt from the manuscript in preparation. 

> In this paper, we introduce a novel method to estimate S by first constructing mirror images Y of the original feature X via Householder reflection. Their difference, W= 1/2 (X-Y), is orthogonal to X and thus T = cov(X, W) is a diagonal matrix.  We define S = \alpha T, and estimate the scalar \alpha by examining the largest eigenvalues of the scaled W.  Such an S=\hat\alpha T leads to much more powerful knockoffs, and its computation requires no optimization.  We call the knockoffs constructed via Householder reflection **reflection knockoff** (ReKo).

The basics of knockoff filter can be found here: https://web.stanford.edu/group/candes/knockoffs/.

## Usage: 
./reko -f filename -c number -t thread -o output_prefix
1. -f input file is a n by p matrix where n is the number samples and p is number of features.
2. -c specify number of copies of reflection knockoff you want to construct.
3. -t specify number of threads. it will use all the avaiable threads by default.
4. -o output prefix.  
