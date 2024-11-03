function [X_sel,H_sel,S,fsamp,L] = data_generator(k)
% EMG data generation
% k = number of sources
% load data
load('SynthMUAP.mat')

%% MUAPs
H = [];
for ii = 1 : length(MUAPs)
    MUAP_ii = MUAPs{1,ii};
    H_ii = [];
    for jj = 1 : length(MUAP_ii)
        H_ii = [H_ii; MUAP_ii{1,jj}];
    end
    H = [H H_ii];
end
%% sFirings
L = 35; % extracted from MUs length
T = SigLength*fsamp;
S = [];
for ii = 1 : k
    locs_ii = sFirings{1,ii};
    s_ii    = zeros(1,T);
    s_ii(locs_ii) = 1;
    s_ii = [zeros(1,L-1) s_ii];
    S_ii = [];
    for ll = 1 : L
        S_ii = [S_ii; s_ii(1, (L-ll+1) : L-ll+T) ];
    end
    S = [S; S_ii];
    
end
%% Data Model X = HS
Hk = H(:,1:k*L);
X  = Hk * S;

%% Select three random channels for testing
idx = [45 59 195];
X_sel = X(idx,:);
H_sel = H(idx,:);
end

