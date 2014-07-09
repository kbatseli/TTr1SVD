%% small demo to illustrate the use of ttr1svd.m and getAtilde.m

clear all;
clc;

% example tensor from "A Constructive Algorithm for Decomposing a Tensor into a Finite Sum of Orthonormal Rank-1 Terms"
A(:,:,1)=[1 4 7 10;2 5 8 11;3 6 9 12];
A(:,:,2)=[ 13 16 19 22; 14 17 20 23; 15 18 21 24];

% store dimensions of A in n
n=size(A);

% compute the TTr1 decomposition of A
[U,S,V,sigmas]=ttr1svd(A);

% check out sigmas
sigmas

% try a rank-2 approximation
[Atilde,s,outerprod,T]=getAtilde(U,S,V,1:2,n);

% orthogonality of outer product factors
outerprod'*outerprod

% compare our approximation error with bound from singular values
[norm(reshape(A-Atilde,[1 prod(n)])) norm(sigmas(3:end))]

% try a different rank-2 approximation
[Atilde2,s2,outerprod2,T2]=getAtilde(U,S,V,[2 4],n);

% compare our approximation error with bound from singular values
[norm(reshape(A-Atilde2,[1 prod(n)])) norm(sigmas([1 3 5 6]))]

% determine nullspace terms
[Atilde3,s3,nullspace,T3]=getAtilde(U,S,V,[5 6],n);

% check whether nullspace tensors are orthogonal to original A
norm(reshape(A,[1 prod(n)])*nullspace)
