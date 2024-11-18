function C = mat2tens_block(C_mat,L_vec,M_vec)
% This function rewrites a  matrix C_mat of size R*sum(Mr)*sum(Lr) rows by K columns into
% a column cell C of R tensors C{r}=Cr (LrxMrxK) 
% L_vec and M_vec holds the R successive values of Lr and Mr, respectively
if length(L_vec)~=length(M_vec)
   error('Input size vectors must have the same length')
end
K=size(C_mat,2);
R=length(L_vec);    
C=cell(R,1);
L_vec2=[0,L_vec];
M_vec2=[0,M_vec];
for r=1:R
 C{r,1}=mat2tens(C_mat(sum(L_vec2(1:r).*M_vec2(1:r))+1:sum(L_vec2(1:r+1).*...
     M_vec2(1:r+1)),:),2,[L_vec2(r+1),M_vec2(r+1),K]);   
    % This transforms LrxMrxK tensor into MrLrxK matrix
end
end