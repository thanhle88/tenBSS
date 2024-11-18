function Mat = kr_part(B,C,partB,partC)
%KR_PART Partition-Wise Kronecker product      
[J M]=size(B);
[K N]=size(C);
if (sum(partB)~=M) 
    error(['Error: a matrix with ',int2str(M),' columns can not be partitioned in such a way'])
end
if (sum(partC)~=N) 
    error(['Error: a matrix with ',int2str(N),' columns can not be partitioned in such a way'])
end
if length(partB)~=length(partC)
     error('Error: the 2 input matrices do not have the same number of blocks')
end

indB=[0 cumsum(partB)];
indC=[0 cumsum(partC)];
indMat=[0 cumsum(partB.*partC)];

 Mat=zeros(J*K,sum(partB.*partC));
 for i=1:length(partC)  
     Mat(:,indMat(i)+1:indMat(i+1))=fast_kron( B(:,indB(i)+1:indB(i+1)) , C(:,indC(i)+1:indC(i+1)));
 end
end
         