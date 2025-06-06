# ReKo
ReKo (Reflection Knockoff) takes an input of feature matrix and produce multiple copies of knockoffs. The basics of knockoff filter can be found here: https://web.stanford.edu/group/candes/knockoffs/. The manuscript detailing ReKo can be found at https://www.biorxiv.org/content/10.1101/2025.01.16.633369v2


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
