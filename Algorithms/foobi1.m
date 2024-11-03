function A = foobi1(X,R,varargin)
%FOOBI1 Fourth-Order-Only Blind Identification of Underdetermined Mixtures
%   A = FOOBI1(X,R) estimates the (JxR) mixing matrix M from J mixtures of
%   R sources, even in underdetermined cases with J<R. The estimate is 
%   obtained by decomposing the observed quadricovariance as in [1]. The 
%   (JxN) observation matrix X holds N samples of J-channel data.
%
%   A = FOOBI1(X,R,options) is used to set the options of the simultaneous
%   diagonalization algorithm cpd3_sgsd:
%
%       options.TolFun = 1e-5       The tolerance used in the simultaneous
%                                   diagonalization step.
%       options.MaxIter = 200       The maximal number of iterations for 
%                                   the simultaneous diagonalization step.
%
%       Other options of cpd3_sgsd can also be set. See the documentation
%       of that function for details.

%   Authors: Michiel Vandecappelle (Michiel.Vandecappelle@kuleuven.be)
%            Lieven De Lathauwer   (Lieven.DeLathauwer@kuleuven.be)
%
%   References:
%   [1] L. De Lathauwer, J. Castaing, J.-F. Cardoso, "Fourth-order cumulant
%   based blind identification of underdetermined mixtures", IEEE
%   Transactions on Signal Processing, Vol. 55, No. 6, Part 2, June 2007, 
%   pp. 2965-2973.
%
%   Version History:
%   - 2021/06/29   MV      Initial Version

% Set options
p = inputParser;
p.addOptional('TolFun', 1e-5);
p.addOptional('MaxIter', 200);
p.KeepUnmatched = true;
p.parse(varargin{:});
fn = [fieldnames(p.Results); fieldnames(p.Unmatched)];
data = [struct2cell(p.Results); struct2cell(p.Unmatched)];
options = cell2struct(data, fn);

% Detect number of mixtures
J = size(X,1);

% Compute quadricovariance matrix
C = cum4(X.');
C = reshape(C,J^2,J^2);

% Compress measurements
[u_C,s_C,~] = svd(C,0);
s_C = diag(s_C);
u_C = u_C(:,1:R);
s_C = s_C(1:R);
H = u_C.*s_C.';

% Compute Pst through rank-1 mapping
P = rank1_mapping(H,[J,J]);

% Check uniqueness condition
if R*(R-1)>J^2*(J-1)^2/2
    scheck = svd(P,0);
    error('foobi1:R', ['R has to satisfy R*(R-1) <= J^2*(J-1)^2/2.'...
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
    F = cpd3_sgsd(M,{F0,F0,C0},options);
    Atilde = u_C.*s_C.'*F{1};
end

% Extract mixing matrix by computing rank-1 approximation for each column 
A = nan(J,R);
for r = 1:R
    [ur,~,~] = svd(reshape(Atilde(:,r),J,J));
    A(:,r) = ur(:,1);
end