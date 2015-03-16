function Atilde=getAtilde(U,sigmas,V,sigmaI,n)
% Atilde=getAtilde(U,sigmas,V,sigmaI,n)
% -------------------------------------
% Reconstructs a tensor as a linear combination of rank terms determined
% from a TTr1 decomposition.
%
% Atilde    =   tensor, the desired tensor from taking a linear combination
%               of rank-1 terms,
%
% U         =   cell, cell of U vectors obtained from the TTr1
%               decomposition,
%
% sigmas    =   vector, contains the singular values of each of the rank-1
%               terms in the order as specified in sigmaI. This is NOT
% 				the entire sigmas vector as returned by ttr1svd.m unless
% 				want a full reconstruction,
%
% V         =   cell, cell of V vectors obtained from the TTr1
%               decomposition,
%
% sigmaI    =   vector, contains indices of rank-1 terms that need to be
%               used in the reconstruction. These are the indices of the
%               desired leaves of the last level in the TTr1-tree, counting from
%               left to right,
%
% n         =   vector, dimensions of the original tensor A.
%
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
Atilde=zeros(prod(n),1);

for j=1:length(sigmaI) % for each rank-1 term    
    tempEnd=endI;
    tempR=r;
    k=sigmaI(j); % the offset with respect to endI for level d-i, initalize to sigmaI

    %% first step is different because we need to multiply with V as well!   
    % determine new offsets
    if mod(k,tempR(end))==0
        k=k/tempR(end);
        l=tempR(end);
    else
        l=mod(k,tempR(end));
        k=ceil(k/tempR(end));
    end
    outerprod=kron(V{tempEnd(end)+k}(:,l),U{tempEnd(end)+k}(:,l));
    
    %% now do remaining steps
    tempEnd=tempEnd(1:end-1);
    tempR=tempR(1:end-1);
    
    for i=1:length(tempEnd)
        % determine new offsets
        if mod(k,tempR(end))==0
            k=k/tempR(end);
            l=tempR(end);
        else
            l=mod(k,tempR(end));
            k=ceil(k/tempR(end));
        end        
        outerprod=kron(outerprod,U{tempEnd(end)+k}(:,l));
        tempEnd=tempEnd(1:end-1);
        tempR=tempR(1:end-1);
    end
    
    %% last step is always first node, so k=1   
    % determine new l-offset
    if mod(k,tempR(end))==0
        l=tempR(end);
    else
        l=mod(k,tempR(end));
    end
    outerprod=kron(outerprod,U{1}(:,l));
    
    Atilde=Atilde+sigmas(j)*outerprod;
end
Atilde=reshape(Atilde,n);

end
