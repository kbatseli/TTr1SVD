Tensor Train Rank-1 decomposition implementation for Matlab&copy;/Octave&copy;
------------------------------------------------------------------------------

The Tensor Train Rank-1 (TTr1) decomposition decomposes an arbitrary tensor A into a unique linear combination of orthonormal rank-1 terms. This allows for a very straightforward determination of a low-rank approximation as well as an easy quantification of the approximation error.

1. Functions
------------

* [U,S,V,sigmas]=ttr1svd(A)

Use this function to compute the TTr1 decomposition.

* [Atilde,sigmas,outerprod,T]=getAtilde(U,S,V,sigmaI,n)

Use this function to compute a rank-R approximation of A from the TTr1 decomposition. Also returns the orthogonal outer product factors.

* demo.m

Small demo that illustrates the use of ttr1svd.m and getAtilde.m.

2. Reference
------------

"A Constructive Algorithm for Decomposing a Tensor into a Finite Sum of Orthonormal Rank-1 Terms"

http://arxiv.org/abs/1407.1593


Authors: Kim Batselier, Haotian Liu, Ngai Wong