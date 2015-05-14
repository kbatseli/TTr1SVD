function ocv=verifyOrtc(A,O)
% ocv=verifyOrtc(A,O)
% -------------------
% This function verifies that each column of O forms a tensor that is
% orthogonal to A.
%
% ocv       =   vector, each ith entry is an inner product of outer product of
%               O{:,i} with A. This should be numerically zero,
%
% A         =   array, d-way array,
%
% O         =   cell, orthogonal tensors computed with orthc.m.
%
% Reference
% ---------
%
% A Constructive Algorithm for Decomposing a Tensor into a Finite Sum of Orthonormal Rank-1 Terms
% http://arxiv.org/abs/1407.1593
%
% 2015, Kim Batselier, Haotian Liu, Ngai Wong

[d,R]=size(O);
ocv=zeros(1,R);

for i=1:R
    temp=O{d,i};
    for k=1:d-1
        temp=kron(temp,O{d-k,i});
    end
    ocv(i)=A(:)'*temp;    
end

end