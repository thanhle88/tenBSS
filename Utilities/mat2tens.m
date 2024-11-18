function [X]=mat2tens(X_mat,mode,size_vec)
%MAT2TENS Tensorization of a matrix (reciprocal of tens2mat)     
if mode<1 || mode >3
    error ('Input argument mode must be a 1 , 2 or 3')
end
I=size_vec(1);
J=size_vec(2);
K=size_vec(3);
if mode==1
    X=permute(reshape(X_mat,K,I,J),[2 3 1]);
elseif mode==2
    X=reshape(X_mat,I,J,K);
elseif mode==3
    X=permute(reshape(X_mat,J,K,I),[3 1 2]);
end
end