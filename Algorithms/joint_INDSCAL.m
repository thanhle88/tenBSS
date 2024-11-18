function [Fac_X,Fac_Y] = joint_INDSCAL(X,Y,R,OPTS)
% Joint INDSCAL decomposition of X (M x M x L) and Y (M2 x M2 x L)
% Author: Thanh Trung Le,
% Vietnam National University, Hanoi
% Email: thanhletrung@vnu.edu.vn

%%
if isfield(OPTS,'rho'),
    rho = OPTS.rho;
else
    rho = 4;
end
rhoA = sqrt(rho); rhoB = rho;
if isfield(OPTS,'Iter_Max'),
    Iter_Max = OPTS.Iter_Max;
else
    Iter_Max = 100;
end
if isfield(OPTS,'Iter_ADMM_Loop'),
    Iter_ADMM_Loop = OPTS.Iter_ADMM_Loop;
else
    Iter_ADMM_Loop = 10;
end

%% 
M = size(X,1); % data dimension

%% Initialization
Factor_Y = cpd(Y,R);

%% Processing

X1 = ten2mat(X,1);
X2 = ten2mat(X,2);
X3 = ten2mat(X,3);
Y1 = ten2mat(Y,1);
Y2 = ten2mat(Y,2);
Y3 = ten2mat(Y,3);

Bk = Factor_Y{1}; 
Bk_old = Bk;
Ek = Bk;
Uk = zeros(M*M,R);
Sk = X3 * pinv(Bk'); 
alpha = (norm(Y(:))/norm(X(:)))^2;

for k = 1 : Iter_Max

 
    %% Stage 1: B_k and K_k
    B_L = Bk;
    B_R = Bk;
    D_l = zeros(M*M,R);
    for l =  1 : Iter_ADMM_Loop
        B_L_old = B_L;
        D_old   = D_l;
        Kk  = Y3 * pinv(khatrirao(B_L,B_R)') ;

        KB_R = khatrirao(Kk,B_R);
        B_L  = (Y1 * KB_R + rho * (Ek - Uk) + rhoB * (B_R - D_l) ) ...
              *(KB_R' * KB_R  + (rho + rhoB)* eye(R) )^(-1);
        
        KB_L = khatrirao(Kk,B_L);
        % B_R  = (Y2 * KB_L  + rhoB * (B_L - D_l)) * (KB_L' * KB_L + rho*eye(R))^(-1);
         B_R = (Y2 * KB_L + rho * (Ek - Uk) + rhoB * (B_L - D_l) ) ...
          *(KB_L' * KB_L  + (rho + rhoB)* eye(R) )^(-1);

        D_l   = D_l + B_L - B_R;
        
        er_t =  cpderr(B_L,B_L_old);
        er_d = norm(D_l-D_old)/norm(D_old); D_old = D_l;
   
        if er_t < 1e-2 && er_d <1e-2 break 
        end
    end
    
    Bk = (B_L+B_R)/2;

    %% AK, E_k, Sigma_k, 

    Ek = (X3' * Sk + rho/alpha*(Uk + Bk) ) * (Sk' *Sk + rho*eye(R))^(-1);   
    Ak  = khatri_rao_inv(Ek,M);
    A_L = Ak;
    A_R = Ak;
    F_l = zeros(M,R);
    for l =  1 : Iter_ADMM_Loop
        A_L_old = A_L;

        Sk  = X3 * pinv(khatrirao(A_L,A_R)') ;

        KA_R = khatrirao(Sk,A_R);
        A_L  = (X1 * KA_R + rhoA * (A_R - F_l) ) ...
               * (KA_R' * KA_R  + (rhoA)* eye(R) )^(-1);
        
        KA_L = khatrirao(Sk,A_L);
        A_R  = (X2 * KA_L  + rhoA * (A_L - F_l)) * (KA_L' * KA_L + rhoA*eye(R))^(-1);
       
        F_l   = F_l + A_L - A_R;
    
        er_t =  cpderr(A_L,A_L_old);
        if er_t < 1e-2 break 
        end
    end
    Ak = (A_L + A_R) /2;
    Ek = khatrirao(A_L,A_R);

    Uk = Uk + Bk - Ek;
    %% Performance 
    er_k = cpderr(Bk,Bk_old); 
    Bk_old = Bk;
    if er_k < 1e-4 break
    end
end
% save
Ak2 = khatri_rao_inv(Bk,M);
Fac_X{1} = Ak2;
Fac_X{2} = Ak;
Fac_X{3} = Sk;
Fac_Y{1} = Bk;
Fac_Y{2} = Bk;
Fac_Y{3} = Kk;

end



