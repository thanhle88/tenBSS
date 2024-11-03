clear;  clc; close all
run_path;

%% load synthetic EMG data 
R = 3;
[EMG,MU,EMG_Source,fsamp,L] = data_generator(R);
% select EMG signals from 2s-5s
T = [fsamp*2:fsamp*5-1];
X_Data = EMG(:,T); 
S = EMG_Source(:,T);
H = MU;  
for r = 1 : R
    H_cell{r} = H(:,[L*(r-1)+1:r*L]);
    Ground_truth{r} = H_cell{r} * S([L*(r-1)+1:r*L],:);
end

%% Algorithms
K = 5;
[~,S_sobi_cp3] = sobi_cp3(X_Data,R,K);
[~,S_ours]     = TCBSS(X_Data,R,[]);

%% PLOTs

fig1 = figure; 
subplot(311); plot(Ground_truth{1,2}(1,:)); 
title('ground truth'); set(gca,'fontsize',15)
subplot(312); plot(S_sobi_cp3(1,:));  
title('sobi-cp3'); set(gca,'fontsize',15) 
subplot(313); plot(-S_ours(1,:));     
title('Proposed type2-BTD'); set(gca,'fontsize',15)
set(fig1, 'units', 'inches', 'position', [0.2 0.5 13 6]);

