function [S,Q]=ttr12tucker(U,sigmas,V,n)
% [S,Q]=ttr12tucker(U,sigmas,V,n)
% -------------------------------
% Converts a tensor A of size n in the TTr1SVD format to the Tucker (HOSVD) format.
% Usually results in a more sparse core S compared to traditional methods
% (e.g. Alternating Least Squares). By setting entries of the sigmas vector
% to zero the Tucker format of a truncated TTr1 will be obtained.
%
% S         =   d-way array, the core tensor in the Tucker decomposition,
%
% Q         =   cell, each Q{k} contains an orthogonal matrix of the
%               kth mode outer product vectors, 
%
% U         =   cell, contains the U vectors of each of the SVDs in the
%               TTr1 tree,
%
% sigmas    =   vector, contains the final singular values in the linear
%               combination of rank-1 terms.
%
% V         =   cell, contains the V vectors of each of the SVDs in the
%               TTr1 tree,
%
% n         =   vector, size of the original tensor A.
%
% Reference
% ---------
%
% A Constructive Algorithm for Decomposing a Tensor into a Finite Sum of Orthonormal Rank-1 Terms
% http://arxiv.org/abs/1407.1593
%
% 2015, Kim Batselier, Haotian Liu, Ngai Wong

d=length(n);
indices=leave2ind(find(sigmas),n); % only handle nonzero sigmas
nonzerosigmas=sigmas(find(sigmas));

Ut=cell(1,d);
% Concatenate all U and V vectors along each mode
for i=d:-1:2
    for j=min(indices(:,2*(i-2)+1)):max(indices(:,2*(i-2)+1))
        I=find(indices(:,2*(i-2)+1)==j);
        Ut{d-i+1}=[Ut{d-i+1} U{j}(:,intersect(indices(I,2*(i-2)+2),indices(I,2*(i-2)+2)))];
        if i==2
            Ut{d}=[Ut{d} V{j}(:,intersect(indices(I,2*(i-2)+2),indices(I,2*(i-2)+2)))];
        end
    end
end

for i=1:d,
    [Q{i},R{i}]=qr(Ut{i},0); % economic QR factorizations of the concatenated U's and V's
    r(i)=rank(Ut{i});
    % distribute the R vectors in a new R cell according to the indices
    % in order for the summation over the nonzero sigmas to be correct
    temp=R{i};
    R{i}=zeros(size(temp,1),length(nonzerosigmas));
    if i==d
        counter=1;
        Rveccounter=1;
        while counter<=size(indices,1)
            I=intersect(find(indices(:,1)==indices(counter,1)),find(indices(:,2)==indices(counter,2)));
            counter=I(end)+1; 
            R{i}(:,I)=temp(:,Rveccounter)*ones(1,length(I));            
            Rveccounter=Rveccounter+1;
        end
    else
        counter=1;
        Rveccounter=1;
        while counter<=size(indices,1)
            I=intersect(find(indices(:,end-(i-1)*2-1)==indices(counter,end-(i-1)*2-1)),find(indices(:,end-(i-1)*2)==indices(counter,end-(i-1)*2)));
            counter=I(end)+1; 
            R{i}(:,I)=temp(:,Rveccounter)*ones(1,length(I));            
            Rveccounter=Rveccounter+1;
        end
    end
end

% construct the core tensor
S=zeros(prod(r),1);
for k=1:length(nonzerosigmas)
    temp=nonzerosigmas(k)*kron(R{d}(:,k),R{d-1}(:,k));
    for i=d-2:-1:1
        temp=kron(temp,R{i}(:,k));
    end
    S=S+temp;    
end
tol=max(n)*eps(max(nonzerosigmas));
S(abs(S)<tol)=0;
S=reshape(S,r);
end
