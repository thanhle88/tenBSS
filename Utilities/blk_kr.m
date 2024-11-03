function C = blk_kr(A,B,L_vec)

C = [];
Nc = cumsum([0,L_vec]);
for r = 1 : length(L_vec)
    C_ii = kron(A(:,Nc(r)+1:Nc(r+1)),B(:,Nc(r)+1:Nc(r+1)));
    C = [C C_ii];
end

end