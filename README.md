Tensor Train Rank-1 Singular Value and Symmetric Eigenvalue Decomposition for Matlab&copy;/Octave&copy;
-------------------------------------------------------------------------------------------------------

The Tensor Train Rank-1 singular value decomposition (TTr1SVD) decomposes an arbitrary tensor A into a unique linear combination of orthonormal rank-1 terms. The Tensor Train Rank-1 symmetric eigenvalue decomposition (TTr1SED) does the same for a tensor that is symmetric in the last 2 modes of equal dimension. This allows for a very straightforward determination of a low-rank approximation as well as an easy quantification of the approximation error.

1. Functions
------------

* [U,S,V,sigmas]=ttr1svd(A) or [U,S,V,sigmas]=ttr1sed(A) 

Use this function to compute the TTr1SVD(/SED) decomposition.

* Atilde=getAtilde(U,sigmas,V,sigmaI,n)

Use this function to compute a rank-R approximation Atilde of A from the TTr1SVD(/SED) decomposition.

* demo.m

Small demo that illustrates the use of ttr1svd.m, getAtilde.m and leave2ind.m.

* O=orthc(A) or O=orthc(A,tol)

This function computes all outer vector products that form tensors orthogonal to A. 

* ocv=verifyOrtc(A,O)

Use this function to verify that all tensors that can be formed with vectors in O are orthogonal to A.

2. Reference
------------

"A Constructive Algorithm for Decomposing a Tensor into a Finite Sum of Orthonormal Rank-1 Terms"

http://arxiv.org/abs/1407.1593


Authors: Kim Batselier, Haotian Liu, Ngai Wong
