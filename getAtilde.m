function [Atilde,sigmas,outerprod,T]=getAtilde(U,S,V,sigmaI,n)
% [Atilde,sigmas,outerprod,T]=getAtilde(U,S,V,sigmaI,n)
% -----------------------------------------------------
% Reconstructs a tensor as a linear combination of rank terms determined
% from a TTr1 decomposition. Also returns the singular values and
% orthogonal outer product factors that make up each of the rank-1 terms.
%
% Atilde    =   tensor, the desired tensor from taking a linear combination
%               of rank-1 terms,
%
% sigmas    =   vector, contains the singular values of each of the rank-1
%               terms in the order as specified in sigmaI,
%
% outerprod =	matrix, each column of this matrix correponds with a
%               unit-norm outer product factor in the order as specified in
%               sigmaI. reshape(outerprod(:,k),n) is the rank-1 tensor
%               corresponding with the kth outer product factor,
%
% T 		=	cell, contains for each rank-1 term the corresponding intermediate
%               singular values, U and V vectors,
%
% polysys   =   cell containing coefficients and monomials exponents of the
%               set of polynomial equations.
%
% U         =   cell, cell of U vectors obtained from the TTr1
%               decomposition,
%
% S         =   cell, cell of singular values obtained from the TTr1
%               decomposition,
%
% V         =   cell, cell of V vectors obtained from the TTr1
%               decomposition,
%
% sigmaI    =   vector, contains indices of rank-1 terms that need to be
%               used in the reconstruction,
%
% n         =   vector, dimensions of the original tensor A.

% Reference
% ---------
%
% A Constructive Algorithm for Decomposing a Tensor into a Finite Sum of Orthonormal Rank-1 Terms
% http://arxiv.org/abs/1407.1593
%
% 2014, Kim Batselier, Haotian Liu, Ngai Wong

% get indices of end node for each level based on size of original tensor n
d=length(n);                            % order of the tensor
r=ones(1,d);                            % nodes are paired in groups of r(i) on level i-1
nodesperlevel=ones(1,d);                % level i-1 contains nodesperlevel(i) nodes
endI=ones(1,d);                         % endI(i) is the index of the last SVD of level i-1
for i=2:d
    r(i)=min(n(i-1),prod(n(i:end)));
    nodesperlevel(i)=prod(r(1:i));
    endI(i)=sum(nodesperlevel(1:i));
end

%last endI points to last SVD of previous level       
endI=endI(1:end-2); 

% initialize output
T=cell(1,length(sigmaI));
outerprod=cell(1,length(sigmaI));
Atilde=zeros(prod(n),1);
sigmas=zeros(1,length(sigmaI));

% k(i) containst the offset of the ith term with respect to endI for
% level d-i, initalize to sigmaI
k=sigmaI;
% l tells us which element of S{k(i)}(:,l(i)) we need to choose
l=zeros(1,length(sigmaI));

%% first step is different because we need to multiply with V as well!
for j=1:length(sigmaI)
    
    % determine offset for each node
    if mod(k(j),r(end))==0
        k(j)=k(j)/r(end);
        l(j)=r(end);
    else
        l(j)=mod(k(j),r(end));
        k(j)=ceil(k(j)/r(end));
    end
    
    T{j}{1}=S{endI(end)+k(j)}(l(j));
    T{j}{2}=V{endI(end)+k(j)}(:,l(j));
    T{j}{3}=U{endI(end)+k(j)}(:,l(j));
    
    sigmas(j)=S{endI(end)+k(j)}(l(j));
    outerprod{j}=kron(V{endI(end)+k(j)}(:,l(j)),U{endI(end)+k(j)}(:,l(j)));
end

%% now do remaining steps
endI=endI(1:end-1);
r=r(1:end-1);
Tcounter=4;
for i=1:length(endI)
    
    for j=1:length(sigmaI)
        % determine offset for each node
        if mod(k(j),r(end))==0
            k(j)=k(j)/r(end);
            l(j)=r(end);
        else
            l(j)=mod(k(j),r(end));
            k(j)=ceil(k(j)/r(end));
        end
        
        T{j}{Tcounter}=S{endI(end)+k(j)}(l(j));
        T{j}{Tcounter+1}=U{endI(end)+k(j)}(:,l(j));
        sigmas(j)=sigmas(j)*S{endI(end)+k(j)}(l(j));
        outerprod{j}=kron(outerprod{j},U{endI(end)+k(j)}(:,l(j)));
    end    
    Tcounter=Tcounter+2;
    endI=endI(1:end-1);
    r=r(1:end-1);

end

%% last step is always first node
for j=1:length(sigmaI)
    
    % determine offset for each node
    if mod(k(j),r(end))==0
        l(j)=r(end);
    else
        l(j)=mod(k(j),r(end));
    end
    
    T{j}{Tcounter}=S{1}(l(j));
    T{j}{Tcounter+1}=U{1}(:,l(j));
    sigmas(j)=sigmas(j)*S{1}(l(j));
    outerprod{j}=kron(outerprod{j},U{1}(:,l(j)));
end

% add R rank-1 terms and reshape into tensor
outerprod=cell2mat(outerprod);
Atilde=reshape(sum(outerprod*diag(sigmas),2),n);

end

