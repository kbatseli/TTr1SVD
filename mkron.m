function y=mkron(varargin)
% y=mkron(x1,x2,x3,...)
% ---------------------
% Returns a Kronecker product of multiple input matrices x1,x2,x3... Needs
% at least 2 input arguments.
%
% y         =   matrix, Kronecker product of x1,x2,x3..
%
% x1,2,3,.. =   matrix, matrix of arbitrary size.

% Reference
% ---------
%
% 2014, Kim Batselier

n=length(varargin);

if n<2
    error('Need at least 2 input arguments.')
end

y=varargin{1};
for i=2:n
    y=kron(y,varargin{i});
end

end