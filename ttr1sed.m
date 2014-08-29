function [U,S,V,sigmas]=ttr1sed(A)
% [U,S,V,sigmas]=ttr1sed(A)
% -------------------------
% Tensor Train rank-1 symmetric eigenvalue decomposition. Decomposes a
% tensor A that is symmetric in the 2 last modes into a linear combination
% of orthonormal rank-1 terms. Returns the orthogonal
% vectors U,V and weights S from each of the SVD/EIGs in the TTr1 tree.
% Use the function getAtilde.m to obtain a low rank approximation using the
% U,S,V obtained from this function.
%
% U         =   cell, contains the U vectors of each of the SVDs in the
%               tree,
%
% S         =   cell, contains the singular values S of each of the SVDs
%               in the tree,
%
% V         =   cell, contains the V vectors of each of the SVDs in the
%               tree,
%
% sigmas    =   vector, contains the final weights in the linear
%               combination of rank-1 terms.
%
% Reference
% ---------
%
% A Constructive Algorithm for Decomposing a Tensor into a Finite Sum of Orthonormal Rank-1 Terms
% http://arxiv.org/abs/1407.1593
%
% 2014, Kim Batselier, Haotian Liu, Ngai Wong


n=size(A);
r=zeros(1,length(n)-1);
for i=1:length(n)-1
    r(i) = min(n(i),prod(n(i+1:end)));
end


totalsvd=1;
svdsperlevel=zeros(1,length(r));   % add length 1 for the first level
svdsperlevel(1)=1;  % first level
for i=2:length(r)
    svdsperlevel(i)=prod(r(1:i-1));
    totalsvd=totalsvd+svdsperlevel(i); 
end
nleaf=prod(r);

U=cell(1,totalsvd);
S=cell(1,totalsvd);
V=cell(1,totalsvd);

[Ut St Vt]=svd(reshape(A,[n(1),prod(n(2:end))]),'econ');
U{1}=Ut;
S{1}=diag(St);
V{1}=Vt;
counter=2; % this counter keeps track during the iterations which V{i} we are computing. This is a linear counter that counts breadth-first
whichvcounter=1;    % this counter keeps track during the iterations of which V we are taking the svd

for i=1:length(r)-1           % outer loop over the levels
    for j=1:prod(r(1:i))      % inner loop over the number of svds for this level 
        if rem(j,r(i)) == 0
            col=r(i);
        else
            col=rem(j,r(i));
        end
        if counter > sum(svdsperlevel(1:end-1))
        % last level generation, use eig now			
			[Vt St]=eig(reshape(V{whichvcounter}(:,col),[n(i+1),prod(n(i+2:end))]));
			Ut=Vt;
        else
			[Ut St Vt]=svd(reshape(V{whichvcounter}(:,col),[n(i+1),prod(n(i+2:end))]),'econ');
		end
        U{counter}=Ut;
        S{counter}=diag(St);
        V{counter}=Vt;
        counter=counter+1;
        if rem(j,length(S{whichvcounter}))==0
            V{whichvcounter}=[];
            whichvcounter =  whichvcounter+1;
        end
    end
%     whichvcounter = whichvcounter+1;
end

Slevel=cell(1,length(r));   % cat each level singular values into 1 vector
counter=1;
for i=1:length(r),
    for j=1:svdsperlevel(i),
        Slevel{i}=[Slevel{i}; S{counter}];
        counter=counter+1;
    end
end

for i=1:length(r),             % make all singular value vectors the same size (number of leaves)
    Slevel{i}=kron(Slevel{i}, ones(nleaf/length(Slevel{i}),1));
end

sigmas=ones(nleaf,1);         % output singular values at each leaf
for i=1:length(r),
    sigmas=sigmas.*Slevel{i};
end

end
