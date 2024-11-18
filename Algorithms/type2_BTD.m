function [A_es,er_A,er_X] = type2_BTD(X,L_vec,OPTS)

% Author: Thanh Trung Le,
% Vietnam National University, Hanoi
% Email: thanhletrung@vnu.edu.vn

if nargin <= 2
    flag = 0;             % without performance estimation part
    er_A = []; er_X = [];
else
    flag = 1;             % with performance estimation part
    X_true = OPTS.X_true;
    A_true = OPTS.A_true;
end

if isfield(OPTS,'data_type')
    data_type = OPTS.data_type;
else
    data_type = 'complex';
end
if isfield(OPTS,'iter_max')
    iter_max = OPTS.iter_max;
else
    iter_max = 500;
end
if isfield(OPTS,'penalty')
    rho = OPTS.penalty;
else
    rho = 1e-3;
end


I = size(X,1);
K = size(X,3);
R = length(L_vec);
N = sum(L_vec);
Nc  = cumsum([0,L_vec]);

%% Processing ...


% Initialize at random
if strcmp(data_type,'complex')
    A_es = randn(I,N) + i*randn(I,N);
    B_es = randn(I,N) + i*randn(I,N);
    D_es = cell(R,1);
    for r = 1 : R
        D_es{r} = randn(L_vec(r),L_vec(r),K) + i*randn(L_vec(r),L_vec(r),K);
    end
else
    A_es = randn(I,N);
    B_es = randn(I,N);
    D_es = cell(R,1);
    for r = 1 : R
        D_es{r} = randn(L_vec(r),L_vec(r),K);
    end
end
U  = zeros(I,N);
Z  = zeros(I,N);

U_old = U;
alpha_old = 1;
A_old = A_es;

X_1 = ten2mat(X,1);
X_2 = ten2mat(X,2);
X_3 = ten2mat(X,3);

er_X = []; er_A = [];

%% Processing ...
for ii =  1 :  iter_max
    % estimate of A
    W_A = cell(1,R);
    for r = 1 : R
        B_r       = B_es(:,Nc(r)+1:Nc(r+1));
        W_A_r_ten = tmprods(D_es{r},B_r,2);
        W_A{r}    = ten2mat(W_A_r_ten,1)';
    end
    W_A    = cell2mat(W_A);
    A_es   = (X_1*W_A + rho*(B_es-U))*pinv(W_A'*W_A + rho*eye(N));

    % estimate of B
    W_B = cell(1,R);
    for r = 1 : R
        A_r = A_es(:,Nc(r)+1:Nc(r+1));
        W_B_r_ten = tmprods(D_es{r},A_r,1);
        W_B{r}    = ten2mat(W_B_r_ten,2)';
    end
    W_B = cell2mat(W_B);
    B_es = (X_2*W_B + rho*(U-A_es))*pinv(W_B'*W_B + rho*eye(N));

    % estimate of G
    AB  = blk_kr(B_es,A_es,L_vec);
    GG  = pinv(AB) * X_3.';
    ER_X = X_3.' - AB * GG;
    D_es = mat2tens_block(GG,L_vec,L_vec);

    U = Z + B_es - A_es;
    alpha = (1 + sqrt(1+4*alpha_old^2))/2;
    Z = U + (alpha_old-1)/alpha*(U-U_old);

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
            X_r = X_r + tmprods(tmprods(D_es{r},A_r,1),B_r,2);
        end
        ER = X_r - X_true;
        er_X_ii = norm(ER(:))/norm(X_true(:));
        er_A_ii = solve_blockperm(A_es,A_true,L_vec);
        er_X = [er_X, er_X_ii];
        er_A = [er_A, er_A_ii];
    end

end

end

