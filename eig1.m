function [eigvec, eigval, eigval_full] = eig1(A, c, isMax, isSym)
%EIG1 此处显示有关此函数的摘要
%   此处显示详细说明
if nargin < 2
    c = size(A,1);
    isMax = 1;
    isSym = 1;
elseif c > size(A,1)
    c = size(A,1);
end;

if nargin < 3
    isMax = 1;
    isSym = 1;
end;

if nargin < 4
    isSym = 1;
end;

if isSym == 1
    A = max(A,A');
end;
try
    [v, d] = eig(A);
    d = diag(d);
    %d = real(d);
catch
    if isMax == 0 
        [v, d] = eigs(sparse(A), c, 'sa', struct('tol', 1e-5'));
    else
        [v, d] = eigs(sparse(A), c, 'la', struct('tol', 1e-5'));
    end
end

if isMax == 0
    [d1, idx] = sort(d);
else
    [d1, idx] = sort(d,'descend');
end;
idx1 = idx(1:c);
eigval = d(idx1);
eigvec = v(:,idx1);

eigval_full = d(idx);

end

