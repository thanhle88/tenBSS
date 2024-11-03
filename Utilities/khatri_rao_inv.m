function A = khatri_rao_inv(B,n)
% recover A from B = khatri_rao(A,A)
r = size(B,2);
A = [];
for ii = 1 : r
    b1 = B(:,ii);
    B1 = reshape(b1,[n n]);
    [u1,sigma1,v1] = svd(B1,'econ');
    u1 = u1(:,1);
    A(:,ii) = u1;
end

end