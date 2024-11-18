function C_mat = tens2mat_block(C)
% This function rewrites a cell C of R tensors C{r}=Cr (LrxMrxK) into a matrix C_mat
% of size R*sum(Mr)*sum(Lr) rows by K columns.
% All tensors Cr must have the same 3rd dimension K.
if min(size(C))~=1
   error('The cell C must have one dimension equal to 1')
end
R=max(size(C));    % The R tensors are stacked in the cell C (either in column of row format)
for r=1:R
 D{r,1}=tens2mat(C{r},2);   % This transforms LrxMrxK tensor into MrLrxK matrix
end
C_mat=cell2mat(D);
end      
       