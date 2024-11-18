function [H_es,S_es] = TCBSS(Data,R,OPTS)
% Author: Thanh Trung Le,
% Vietnam National University, Hanoi
% Email: thanhletrung@vnu.edu.vn

if isfield(OPTS,'tau'), % number of time lags
    tau = OPTS.tau;
else
    tau = 5;
end
if isfield(OPTS,'MU_length'),
    L = OPTS.MU_length;
else
    L = 35;
end

X = Data;
M = size(X,1);
T = size(X,2);

L_      = 10; % extra observations
X_block = [];
for ii = L_+1 : T
    x_ii = [];
    for m = 1 : M
        x_m = X(m, ii-L_+1:ii);
        x_ii = [x_ii,x_m];
    end
    X_block = [X_block x_ii'];
end
%
RX    = corrX(X,tau);
RX2   = corrX(X_block,tau);

L_vec  = 3*ones(1,R);

%% Tensor decomposition
[H_es,~,~] = BTD2_IALM(RX,L_vec);
S_es = pinv(H_es)*X;

end


function Rx = corrX(X,tau)
[M,P] = size(X);
% Make mixtures zero mean
for m = 1 : M
    X(m,:) = X(m,:)- mean(X(m,:));
end

if exist('tau')~=1 || isempty(tau);tau=2*M;end
if length(tau)==1; tau=1:tau; end
K = length(tau);

Rx = zeros(M,M,K);
for t = 1:K
    Rx(:,:,t) = X(:,tau(t)+1:end) * X(:,1:end-tau(t))' / (P-tau(t));
end
% Symmetrization of covariance matrices
Rx = 0.5*(Rx+conj(permute(Rx,[2 1 3])));
end


function [A_es,er_A,er_X] = BTD2_IALM(X,L_vec,OPTS)

if nargin <= 2
    flag = 0;             % without performance estimation part
    er_A = []; er_X = [];
else
    flag = 1;      % with performance estimation part
    X_true = OPTS.X_true;
    A_true = OPTS.A_true;
end

I = size(X,1);
K = size(X,3);

R  = length(L_vec);
N  = sum(L_vec);
Nc = cumsum([0,L_vec]);
N2c = cumsum([0,L_vec.^2]);

%% Processing ...
Iter_Max = 100;
rho = 1;
% Initialize
A_es = randn(I,N);
B_es = randn(I,N);
U    = zeros(I,N);
Z    = zeros(I,N);

U_old = U;
Z_old = Z;
alpha_old = 1;
A_old = A_es;

D_es = cell(R,1);
for r = 1 : R
    D_es{r} =  randn(L_vec(r),L_vec(r),K);
end

X_1 = ten2mat(X,1);
X_2 = ten2mat(X,2);
X_3 = ten2mat(X,3);

er_X = []; er_A = [];

%% Processing ...
for ii =  1 : Iter_Max
    % estimate of A
    W_A = [];
    for r = 1 : R
        B_r = B_es(:,Nc(r)+1:Nc(r+1));
        W_A_r_ten = tmprod(D_es{r},B_r,2);
        W_A_r = ten2mat(W_A_r_ten,1)';
        W_A = [W_A W_A_r];
    end
    W_A = W_A';
    A_es = (X_1*W_A' + rho*(B_es-U) ) * pinv(W_A*W_A' + rho*eye(N));
    A_cell = mat2cell(A_es,size(A_es,1),L_vec);
    for r=1:R
        [A_cell{r},~] = qr(A_cell{r},0);
    end
    A_es = cell2mat(A_cell);

    % estimate of B
    W_B = [];
    for r = 1 : R
        A_r = A_es(:,Nc(r)+1:Nc(r+1));
        W_B_r_ten = tmprod(D_es{r},A_r,1);
        W_B_r = ten2mat(W_B_r_ten,2)';
        W_B = [W_B W_B_r];
    end
    W_B = W_B';
    B_es   = (X_2*W_B' + rho*(U-A_es) ) * pinv(W_B*W_B' + rho*eye(N));
    B_cell = mat2cell(B_es,size(B_es,1),L_vec);
    for r=1:R
        [B_cell{r},~] = qr(B_cell{r},0);
    end
    B_es = cell2mat(B_cell);

    % estimate of G
    AB = blk_kr(B_es,A_es,L_vec);
    GG = X_3 * pinv(AB');
    for r = 1 : R
        G_r  = GG(:,N2c(r)+1:N2c(r+1));
        D_es{r} = mat2ten(G_r,[L_vec(r),L_vec(r),K],3);
    end

    U = Z + B_es - A_es;
    alpha = (1 + sqrt(1+4*alpha_old^2))/2;
    Z = U + (alpha_old-1)/alpha*(U-U_old);

    Z_old = Z;
    U_old = U;
    alpha_old = alpha;

    % stop iteration
    if  mod(ii,5) == 0
        er = solve_blockperm(A_es,A_old,L_vec);
        A_old = A_es;
        if er <1e-3 && ii>20
            break;
        end
    end

    %% Performance Evaluation Part
    if flag == 1
        X_r = zeros(I,I,K);
        for r = 1 : R
            A_r = A_es(:,Nc(r)+1:Nc(r+1));
            B_r = B_es(:,Nc(r)+1:Nc(r+1));
            X_r = X_r + tmprod(tmprod(D_es{r},A_r,1),B_r,2);
        end
        ER = X_r - X_true;
        er_X_ii = norm(ER(:))/norm(X_true(:));
        er_A_ii = solve_blockperm(A_es,A_true,L_vec);
        er_X = [er_X, er_X_ii];
        er_A = [er_A, er_A_ii];
    end

end

end

