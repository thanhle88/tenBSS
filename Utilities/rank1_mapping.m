function P = rank1_mapping(E,sz,gramian)
%RANK1_MAPPING Matricized rank-1 detection tensor.
%   P = RANK1_MAPPING(E,SZ) constructs a matricization of the rank-1 detection tensor of size
%   I-by-I-by-J-by-J-by-R-by-R where [I,J] = SZ. This matricization takes symmetries in the last two
%   modes into account by taking only the lower triangular part, which results in a matrix
%   representation P of size I^2*J^2-by-(R+1)*R/2. Each column of P is the matricized result of
%   rank-1 detection mapping phi [1, Theorem 2.1] applied to each combination of columns of the
%   I*J-by-R matrix E, i.e., P is a matrix defined as [1,2]:
%
%       P = [phi(Es,Et)]    with 1<=t<=s<=R
%
%   Es and Et denote the matricized s-th and t-th column of E and have dimensions IxJ.
%
%   G = RANK1_MAPPING(E,SZ,TRUE) computes the Gramian of the rank-1 detection device, i.e., G =
%   0.25*P'*P, which has dimensions (R+1)*R/2-by-(R+1)*R/2 and can be computed more efficiently. 
    
%   Authors:    Mikael Sorensen         
%               Nico Vervliet           (Nico.Vervliet@kuleuven.be)
%               Martijn Boussé          (Martijn.Bousse@kuleuven.be)
%               Laurent Sorber          
%               Lieven De Lathauwer     (Lieven.DeLathauwer@kuleuven.be) 
%
%   References:
%   [1] L. De Lathauwer, "A Link between the Canonical Decomposition in Multilinear Algebra and
%       Simultaneous Matrix Diagonalization," SIAM J. Matrix Anal. Appl., Vol. 28, No. 3, 2006,
%       pp. 642-666.
%   [2] I. Domanov, L. De Lathauwer, "Canonical Polyadic Decomposition of Third-Order Tensors:
%       Reduction to Generalized Eigenvalue Decomposition," SIAM J. Matrix Anal. Appl., Vol. 35,
%       No. 2, 2014, pp.  636-660.
%        
%   Version history:
%   - 2022/01/20    NV      Added Gramian computation for non-othogonal E.    
%   - 2022/01/19    NV      Simplified code and added Gramian. 
%   - 23/06/2020    MB      Updated version.
%   - 23/06/2020    MS      Original version.

    % Check inputs.       
    if length(sz) ~= 2
        error('rank1_mapping:sz','length(sz) should be equal to 2.')
    end
    if size(E,1) ~= prod(sz)
        error('rank1_mapping:E','size(E,1) should be equal to prod(sz).')
    end
    if nargin < 3
        gramian = false;
    end 
        
    R = size(E,2);
    if ~gramian 
        % Construct rank-1 detection device matrix for each column.
        P = zeros(prod(sz)^2,R*(R+1)/2);
        E = reshape(E, [sz,R]);
        [s,t] = ind2sub([R,R],find(tril(ones(R))));
        for i = 1:numel(s)
            P(:,i) = phi(E(:,:,s(i)),E(:,:,t(i)));
        end 
    else
        EHE = E'*E;
        isorth = norm(EHE-eye(R),'fro')/norm(E,'fro') < numel(EHE)*2*eps;
        
        % Compute the Gramian 0.25*P'*P for the rank-one detecting device P. The Gramian G is
        % computed as a fourth-order tensor in which an element is given by 
        % 
        %     G(r,s,t,u) = <Er,Et>*<Es,Eu> + <Er,Eu>*<Es,Et> - <Eu'*Es,Er'*Et> - <Et'*Es,Er'*Eu>, 
        %
        % where En denotes the matricized nth column of E.

        % First, create a mask for extracting the lower triangular part of a symmetric matrix.
        mask = tril(true(R));        
        mask = mask(:);
        
        % Compute terms one and two.
        if isorth 
            P12 = ones(R*(R+1)/2,1);
            P12(mask) = 2;
            P12 = diag(P12);
        else 
            W = E'*E;
            P12 = reshape(W,[R,1,R,1]) .* reshape(W,[1,R,1,R]);
            P12 = P12 + reshape(W,[R,1,1,R]) .* reshape(W,[1,R,R,1]);   
            P12 = reshape(P12,R*R,R*R);
            P12 = P12(mask,mask);            
        end         
        % Compute terms three and four.
        Et = reshape(E,sz(1),sz(2),R);
        E1 = reshape(E,sz(1),[]);
        P34 = zeros(sz(2),sz(2),R,R);
        for r = 1:R
            X = Et(:,:,r)'*E1;
            P34(:,:,:,r) = reshape(X,sz(2),sz(2),R);
        end
        P34 = reshape(P34, [sz(2)^2, R^2]);
        P34 = reshape(P34'*P34,[R R R R]);
        P34 = -reshape(permute(P34,[1 4 2 3]),[R^2 R^2]) ...
              -reshape(permute(P34,[1 4 3 2]),[R^2 R^2]);
        % Create final result.
        P = P12 + P34(mask,mask);
    end 
end

function S = phi(X,Y,reduced)
%PHI3 Rank-1 detection mapping.
%   S = PHI(X,Y) returns a vector representation of the rank-1 detection mapping phi(X,Y) of size
%   I-by-I-by-J-by-J defined element-wise as: phi(X,Y)(i,j,k,l) = xij ykl + yij xkl - xkj yil - ykj
%   xil, for matrices X and Y of size I-by-J, see [1, Theorem 2.1].
%
%   S = PHI(X,Y,true) returns a vector representation of the unique elements of the rank-1 detection
%   mapping phi(X,Y).

%   Authors:    Mikael Sorensen         
%               Nico Vervliet           (Nico.Vervliet@kuleuven.be)
%               Martijn Boussé          (Martijn.Bousse@kuleuven.be)
%               Lieven De Lathauwer     (Lieven.DeLathauwer@kuleuven.be) 
%
%   References:
%   [1] L. De Lathauwer, "A Link between the Canonical Decomposition in Multilinear Algebra and
%       Simultaneous Matrix Diagonalization," SIAM J.  Matrix Anal. Appl., Vol. 28, No. 3, 2006,
%       pp. 642-666.
%
%   Version history
%   - 2022/01/19    NV      Simplified code.
%   - 2020/06/23    MS      Original version.
%   - 2020/06/23    MB      Updated version.

    if any(size(X) ~= size(Y))
        error('phi:X','size(X) should be equal to size(Y).')
    end
    X2 = permute(X,[3 4 1 2]);
    Y2 = permute(Y,[3 4 1 2]);

    P1 = X.*Y2 + X2.*Y;
    P2 = permute(P1,[3 2 1 4]);
    D = P1 - P2;
    
    if nargin > 2 && reduced
        % Select unique elements
        [I,J] = size(X);
        TJ = tril(true(J, J), -1);
        TI = tril(true(I, I), -1);
        TIJ = bsxfun(@and, reshape(TJ, [1 J 1 J]), reshape(TI, [I 1 I 1]));
        TIJ = TIJ(:);
        S = D(TIJ);
    else
        S = D(:);
    end

end