% demo script to test bcdLM_alsls algorithm.
% The purpose of this demo is to show how to exploit several initializations.
% The function bcdLM_alsls has 2 kinds of optional input arguments to deal with initialization:
% - Ninit : the algorithm will generate Ninit different starting points and for
%   each starting point, it will stop when the stop criterion is satisfied 
%  (input parameters Tol1 and MaxIt1 to select tolerance and max number of iterations).
%  Output parameters are the ones associated to the best starting point.
% - A,B,C: if one or more matrices are provided, the latter are used as a starting point
%   and this is the only starting point that will be tried, no matter the value of Ninit
%  (it is for instance useful when one already knows estimates close to the solution).

% When one uses e.g. Ninit=10, it can take time before the algorithm stops because
% each starting point will be tried until the convergence criterion is satisfied.
% Thus, it can be useful in practice to try Ninit=10 starting points for a limited
% number of iterations, say MaxItInit (less than MaxIt1) and then go on until convergence
% with only the best initialization obtained among the Ninit tried.
% This demo show how to to this.

clear all
close all
clc

%**********************************************
%--- Choose PARAMETERS of the DEMO
%**********************************************
    %--- Data parameters ----
    data_type='complex';       % choose 'real' or 'complex' to select the kind of data to generate
    I=14;                      % Dimensions  of the tensor
    J=9;
    K=30;
    R=5;                       % number of components
    L=3;                       % rank of the R component matrices Ar
    M=2;                       % rank of the R component matrices Br
    L_vec=L*ones(1,R);
    M_vec=M*ones(1,R);
    SNR=40;                    % SNR [dB], choose SNR=inf for a noise-free model
    power_max=3;               % Each of the R terms is normalized and then multiplied by a value between 1 and power_max. 
                               % The highest this value, the more difficult the problem is (swamps are more likely)
    power_vec=linspace(1,power_max,R);   % Holds the set of R amplitude values
       

    
%***************************************************
%---- Build Loading matrices and observed tensor-----
%**************************************************** 
A_cell=cell(1,R);
B_cell=cell(1,R);
C_cell=cell(1,R);
if strcmp(data_type,'real')==1
    for r=1:R
    A_cell{r}=randn(I,L);B_cell{r}=randn(J,M);C_cell{r}=randn(L,M,K);
    end
elseif strcmp(data_type,'complex')==1
    for r=1:R
    A_cell{r}=randn(I,L)+j*randn(I,L);B_cell{r}=randn(J,M)+j*randn(J,M);C_cell{r}=randn(L,M,K)+j*randn(L,M,K);
    end
end
% Create observed tensor that follows the BCD-(L,M,.) model
Xr=cell(1,R);
X=zeros(I,J,K);
for r=1:R
    XX=zeros(I,J,K);
    for k=1:K
        XX(:,:,k)=A_cell{r}*C_cell{r}(:,:,k)*B_cell{r}.';
    end
    Xr{r}=power_vec(r)*XX/norm(reshape(XX,I*J,K),'fro');   %--> normalize and weight the rth tensor
    X=X+Xr{r};   % --> Build the observed tensor
end
% Add noise
if strcmp(data_type,'real')==1
   Noise_tens=randn(I,J,K);
else
   Noise_tens=randn(I,J,K)+j*randn(I,J,K);
end
sigma=(10^(-SNR/20))*(norm(reshape(X,J*I,K),'fro')/norm(reshape(Noise_tens,J*I,K),'fro'));
X=X+sigma*Noise_tens;


%-------------------------------------------------------------------------------
% Compute the BCD-(L,M,.) with a limited number of iterations for each starting point
% and then use the best starting point to go on until convergence
%-------------------------------------------------------------------------------
Ninit=10;
MaxItInit=10;
% Perform at most MaxItInit iterations for each starting point
[A_init,B_init,C_init,phi,itInit]=bcdLM_alsls(X,R,L,M,[],[],[],MaxItInit,[],[],Ninit);
% Then go on with the best starting point
[A_est,B_est,C_est,phi,it1]=bcdLM_alsls(X,R,L,M,[],[],[],[],[],[],[],A_init,B_init,C_init);
% compute relative errors
[err_A]=solve_blockperm(cell2mat(A_est),cell2mat(A_cell),L_vec);
[err_B]=solve_blockperm(cell2mat(B_est),cell2mat(B_cell),M_vec);
disp(['err on A:         ',num2str(err_A)])
disp(['err on B:         ',num2str(err_B)])
disp(['phi:              ',num2str(phi)])
disp(['itInit:           ',num2str(itInit)])
disp(['it1:              ',num2str(it1)])
