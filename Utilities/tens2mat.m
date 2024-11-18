function [X_mat]=tens2mat(X,mode)
%TENS2MAT Matrix unfoldings of a 3rd order tensor X along a given mode
% INPUTS: - X : tensor of size (IxJxK)
%         - mode = 1 or 2 or 3
% OUTPUTS: X_mat is the matrix unfolding representation of X
% if mode==1:  X_mat is IKxJ  (i=1,...,I, is the slowly varying index)
% if mode==2:  X_mat is JIxK  (j=1,...,J, is the slowly varying index)
% if mode==3   X_mat is KJxI  (k=1,...,K  is the slowly varying index)
[I,J,K]=size(X);
if mode==1
    X_mat=reshape(permute(X,[3 1 2]),I*K,J);
elseif mode==2
    X_mat=reshape(X,J*I,K);
elseif mode==3
    X_mat=reshape(permute(X,[2 3 1]),K*J,I);
else
    error('Input argument mode must be 1, 2 or 3');
end
end
