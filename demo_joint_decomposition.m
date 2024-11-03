clear;close all;clc;
run_path;

%% Experiment Setup
disp('Setup ...')
% Parameters
M = 5;             % No. Sensors
R = 3;             % No. Sources
L = 100;          % No. Tensor Slices
SNR = 10;          % Signal to noise ratio in dB

% Mixing matrix A 
A = randn(M,R);
B = khatrirao(A,A);

X_true = [];
Y_true = [];
for l = 1 : L
    c_l = (randn(R,1)); C(l,:) = c_l'; % variance
    d_l = (randn(R,1)); D(l,:) = d_l'; % kurtosis

    X_true(:,:,l) =  A * diag(c_l) * A';
    Y_true(:,:,l) =  B * diag(d_l) * B';
end
clear c_l d_l l

X  = noisy(X_true,SNR);
Y  = noisy(Y_true,SNR);


%% ALgorithms
disp('Processing ...')

Factor_X_CP = cpd(X,R);
Factor_Y_CP = cpd(Y,R);
OPTS = [];
[Fac_X,Fac_Y] = joint_INDSCAL(X,Y,R,OPTS);

disp('Results ...')

A_es = Fac_X{1}; er_A = cpderr(A,A_es);
B_es = Fac_Y{1}; er_B = cpderr(B,B_es);
A_es_CP = Factor_X_CP{1}; er_A_CP = cpderr(A,A_es_CP);
B_es_CP = Factor_Y_CP{1}; er_B_CP = cpderr(B,B_es_CP);
A_es_B = khatri_rao_inv(B_es_CP,M); er_A_B = cpderr(A,A_es_B);
disp(' ')
fprintf('+ Individual Tensor Decomposition: \n    Tensor 1: Error in A = %f \n    Tensor 2: Error in A = %f  | B = %f \n\n',er_A_CP,er_A_B,er_B_CP);
fprintf('+ Proposed Joint Decomposition:\n   Error in A = %f | B = %f \n\n',er_A,er_B);


