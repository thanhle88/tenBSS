function M = symmetric_null(P, R)
%SYMMETRIC_NULL Compute symmetric matrices in the null space.
%   M = SYMMETRIC_NULL(P,R) computes R symmetric R-by-R matrices in the R dimensional kernel of the
%   tall matrix P which has R*(R+1)/2 columns. The null space of P should have dimension R. Each
%   matrix is stored as frontal slice of the R-by-R-by-R tensor M, and obeys
%
%       mask = tril(true(R,R));   % lower triangular part
%       ind = linspace(1,R*R,R);  % indices of diagonal entries
%       M = reshape(M,R*R,R);     
%       M(ind,:) = M(ind,:)/2;
%       P*M(mask,:)               % should be zero, as M(mask,:) is in the null space. 
%   
%   See also: rank1_mapping

%   Authors:    Nico Vervliet           (Nico.Vervliet@kuleuven.be)
%               Martijn Boussé          (Martijn.Bousse@kuleuven.be)
%               Lieven De Lathauwer     (Lieven.DeLathauwer@kuleuven.be) 
%
%   Version history:
%   - 2022/01/19    NV      Simplified symmetrization.
%   - 2020/07/31    MB+NV   Original version.

    % Check inputs
    if size(P,2) ~= R*(R+1)/2
        error('symmetric_null:P','size(P,2) should equal R*(R+1)/2')
    end

    % Determine null space of P
    [~,~,Vp] = svd(P,0);
        
    % Create lower triangular part of M.
    M = zeros(R^2,R);
    mask = tril(true(R,R));
    M(mask,:) = Vp(:,end-R+1:end);
    % Symmetrize the result to get the full M.
    M = reshape(M, [R R R]);
    M = M + permute(M, [2 1 3]);    
end