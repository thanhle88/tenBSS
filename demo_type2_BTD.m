clc; clear; close all;
 
run_path;

%*******************************************
R = 2;
L_vec = 5*ones(1,R);         % size of the blocks
I = 20;                      % number of rows of the matrix A
K = 100;                     % number of matrices to block-diagonalize
sigma = 1e-1;                % noise level
N=sum(L_vec);
N2c = cumsum([0,L_vec.^2]);
Nc  = cumsum([0,L_vec]);

%*********************************
%% Generate Data
D = zeros(N,N,K);
A = randn(I,N) + i*randn(I,N);
for r=1:R
    D(Nc(r)+1:Nc(r+1),Nc(r)+1:Nc(r+1),:)=randn(L_vec(r),L_vec(r),K);
end
X_true = tmprods(tmprod(D,A,1),A,2);
Noise_tens = randn(I,I,K) + i*randn(I,I,K);
X = X_true + sigma*Noise_tens;

OPTS.X_true = X_true;
OPTS.A_true = A;



%% State-of-the-art 

%% ALS
Tol1   = 1e-6;        % Tolerance
MaxIt1 = 1000;        % Max number of iterations
Tol2   = 1e-4;        % tolerance in refinement stage (after decompression)
MaxIt2 = 50;          % Max number of iterations in refinement stage
Ninit  = 5;           % Number of initializations

%% 
disp('+ ALS')
[A_init,B_init,C_init]=bcdLM_init(X,L_vec,L_vec);  
[A_est,~,~,~,~,~,phi_als]=bcdLM_alsls(X,R,L_vec(1),L_vec(1),'none','on',Tol1,MaxIt1,Tol2,MaxIt2,Ninit,A_init,B_init,C_init);
err_ALS = solve_blockperm(cell2mat(A_est),A,L_vec);
fprintf('       error = %f \n',err_ALS)

disp('+ ALS-LSH')
[A_est_Harshman,~,~,~,~,~,phi_lsh]=bcdLM_alsls(X,R,L_vec(1),L_vec(1),'lsh','on',Tol1,MaxIt1,Tol2,MaxIt2,Ninit,A_init,B_init,C_init);
err_ALS_Harshman = solve_blockperm(cell2mat(A_est_Harshman),A,L_vec);
fprintf('       error = %f \n',err_ALS_Harshman)

disp('+ ALS-ELSC')
[A_est_ELSC,~,~,~,~,~,phi_elsc]=bcdLM_alsls(X,R,L_vec(1),L_vec(1),'elsc','on',Tol1,MaxIt1,Tol2,MaxIt2,Ninit,A_init,B_init,C_init);
err_ALS_ELSC = solve_blockperm(cell2mat(A_est_ELSC),A,L_vec);
fprintf('       error = %f \n',err_ALS_ELSC)


disp('+ Proposed type2-BTD')
OPTS.data_type = 'complex';
[A_es,er_A,er_X] = type2_BTD(X,L_vec,OPTS);
err_our = er_A(end);
fprintf('       error = %f \n',err_our)


%% 
figure; 
semilogy(er_X,'r','LineWidth',2);  
hold on; semilogy(phi_als/norm(X_true(:)),'b','LineWidth',2);
hold on; semilogy(phi_lsh/norm(X_true(:)),'g','LineWidth',2);
hold on; semilogy(phi_elsc/norm(X_true(:)),'k','LineWidth',2);
legend({'proposed','ALS','LSH','ELSC'},'Location','best');

axis([1 100 1e-3 1])
ylabel('Error');
xlabel('Iteration')
title('convergence rate') 
set(gca,'FontSize',20)
