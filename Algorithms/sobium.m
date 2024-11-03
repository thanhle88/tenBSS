function A = sobium(X,K,R,varargin)
%SOBIUM Blind Identification by Simultaneous Matrix Diagonalization
%   A = SOBIUM(X,K) computes the (JxR) mixing matrix A from J mixtures
%   of R sources, even in underdetermined cases with J<R. The (JxN) 
%   observation matrix X holds N J-dimensional measurements. K is a vector 
%   that holds the lag values used for the sample covariance matrices. K 
%   can also be given as a scalar, in which case the lags 0 to K-1 are used. 
%   The number of sources R is estimated by visual inspection and is 
%   assumed to be smaller than or equal to the number of lags. This
%   corresponds to the specific algorithm in [1, Table II].
%
%   A = SOBIUM(X,K,R) computes the mixing matrix with the specified
%   number of sources R. If the method is not explicity provided (see 
%   below), then the approach from Table I in [1] is used if R>K and the
%   approach from Table II is used otherwise.
%
%   A = SOBIUM(X,K,R,options) or A = SOBIUM(X,K,[],options) is used to set 
%   the following option:
%
%       options.Method          The method to be used. The cases 1 and 2
%       [{'auto'},1,2]          correspond to the approaches from Tables I 
%                               and II from [1], respectively. The default 
%                               'auto' uses method 1 if R>K and method 2
%                               if R<=K or R is set to [].

%   Authors: Michiel Vandecappelle (Michiel.Vandecappelle@kuleuven.be)
%            Lieven De Lathauwer   (Lieven.DeLathauwer@kuleuven.be)
%
%   References:
%   [1] L. De Lathauwer, J. Castaing, "Blind identification of 
%   underdetermined mixtures by simultaneous matrix diagonalization", IEEE
%   Transactions on Signal Processing, Vol. 56, No. 3, March 2008, 
%   pp. 1096-1105.
%
%   Version History:
%   - 2021/06/16   MV      Initial Version

% No R given. Use Method 2.
if nargin<3 || isempty(R)
    R = 0;
end

% If K is a scalar, choose 0 to K-1 as lags
if isscalar(K)
    shift = 0:K-1;
else
    shift = K;
    K = length(shift);
end

p = inputParser;
p.addOptional('Method', 'auto');
p.KeepUnmatched = true;
p.parse(varargin{:});
fn = [fieldnames(p.Results); fieldnames(p.Unmatched)];
data = [struct2cell(p.Results); struct2cell(p.Unmatched)];
options = cell2struct(data, fn);

if ischar(options.Method) && strcmpi(options.Method, 'auto')
    if R>K
        options.Method = 1;
    else
        options.Method = 2;
    end
end

if ~isscalar(options.Method) || (options.Method~=1 && options.Method~=2)
    error('sobium:Method','options.Method should either be 1, 2 or "auto".');
end

if R>K && options.Method == 2
    error('sobium:R', 'For Method 2, R should not be greater than K')	
end

% Get J and N from X
J = size(X,1);

% Define shifts and compute covariance tensor

C = scov(X.',shift);

if options.Method == 1 % General case
    % Compute CPD using optimization
    U = cpd(C,R);
    
    % Extract mixing matrix    
    A = U{1};
    
else % Specific case
    
    % Reshape covariance tensor
    C = tens2mat(C,[],3);
        
    % Compress measurements
    [u_C,s_C,~] = svd(C,0);
    s_C = diag(s_C);
    
    if R == 0 % Try to estimate R by inspecting singular values
       semilogy(s_C,'ro');
        R = input('Estimate R by inspecting the singular values:');
        if ~isscalar(R) || R<1 || R>K
           error('sobium:R', 'R has to be a positive integer smaller than K')	
        end
    end
    u_C = u_C(:,1:R);
    s_C = s_C(1:R);
    H = u_C.*s_C.';
    
    % Compute Pst through rank-1 mapping
    P = rank1_mapping(H,[J,J]);
      
    % Check uniqueness condition
    if R*(R-1)>J^2*(J-1)^2/2
        scheck = svd(P,0);
        error('sobium:R', ['R has to satisfy R*(R-1) <= J^2*(J-1)^2/2.'...
         ' The smallest singular value that should be nonzero is equal to %d.'],scheck(end-R));
    end
    
    % Construct auxiliary tensor
    M = symmetric_null(P,R);
    
    % Simultaneous diagonalization
    if R == 1
        Atilde = u_C.*s_C.'*M;
    else
        [F0,~] = eig(M(:,:,1),M(:,:,2));
        F0 = pinv(F0).';
        C0 = reshape(M,[R*R R]).'*conj(kr(F0,F0))/conj((F0'*F0).^2);
        F  = cpd3_sgsd(M,{F0,F0,C0});
        Atilde = u_C.*s_C.'*F{1};
    end

    % Extract mixing matrix by computing rank-1 approximation for each column 
    A = nan(J,R);
    for r = 1:R
        [ur,~,~] = svd(reshape(Atilde(:,r),J,J));
        A(:,r) = ur(:,1);
    end
end

end