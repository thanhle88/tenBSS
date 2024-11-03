%% Setup
clc; clear; close all;
run_path;

% Parameters
J = 3;           % No. Sensors
R = 3;           % No. Sources
N = 1e5;         % No. Samples
SNR = linspace(0,20,5);     

%% Generate sources and mixture matrix
ER_TenSOFO = zeros(length(SNR),1);
ER_SOBIUM  = zeros(length(SNR),1);
ER_FOOBI   = zeros(length(SNR),1);
% Generate data matrix  
filterlength = 20;
f = randn(N-filterlength+1,R); 
g = randn(filterlength,R);
for r = 1:R
    s(:,r) = conv(f(:,r),g(:,r));  
end
S      = (s');
M      = randn(J,R);  
X_true = M*S;         

% Set number of lags
K = 5;
%% Tensor-based BSS Algorithms 
for ii = 1 : length(SNR)

    fprintf(' SNR = %ddB \n', SNR(ii))
    snr_ii = SNR(ii);
    Noise = randn(size(X_true));
    X = X_true + 10^(-snr_ii/20)/norm(Noise,'fro')*norm(X_true,'fro')*Noise; 
    
    %%  FOOBI
    M_FOOBI = foobi1(X,R);
    ER_FOOBI(ii) =  cpderr(M,M_FOOBI);
    %%  SOBIUM
    M_SOBIUM = sobium(X,K,R);
    ER_SOBIUM(ii) = cpderr(M,M_SOBIUM);
    %% TenSOFO
    OPTS.lags = K;
    [M_TENSOFO,Lk]  = tenSOFO(X,R,OPTS);
    ER_TenSOFO(ii)  = cpderr(M,M_TENSOFO);

end

% % Recover the source signals
% S_es_tensofo = pinv(M_TENSOFO)*X;
% S_FOOBI      = pinv(M_FOOBI)*X;
% S_SOBIUM     = pinv(M_SOBIUM)*X;
 
%% Plot 
makerSize    = 14;
numbMarkers  = 50;
LineWidth    = 2;
set(0, 'defaultTextInterpreter', 'latex');
color   = get(groot,'DefaultAxesColorOrder');
red_o   = [1,0,0];
blue_o  = [0, 0, 1];
gree_o  = 'g';%[0, 0.5, 0];
black_o = [0.25, 0.25, 0.25];
blue_n  = color(1,:);
oran_n  = color(2,:);
yell_n  = color(3,:);
viol_n  = color(4,:);
gree_n  = color(5,:);
lblu_n  = color(6,:);
brow_n  = color(7,:);
lbrow_n = [0.5350    0.580    0.2840];
 
fig2 = figure;
k = 10;
K = length(SNR);
KL = 200;

hold on;
d2 = semilogy(1:K,ER_SOBIUM,...
    'linestyle','-','color',gree_o,'LineWidth',LineWidth);
d21 = plot(1:K,ER_SOBIUM,...
 'marker','+','markersize',makerSize,...
   'linestyle','none','color',gree_o,'LineWidth',LineWidth);
d22 = semilogy(1:1,ER_SOBIUM(1),...
    'marker','+','markersize',makerSize,...
    'linestyle','-','color',gree_o,'LineWidth',LineWidth);
 
d4 = semilogy(1:K,ER_FOOBI,...
    'linestyle','-','color',black_o,'LineWidth',LineWidth);
d41 = plot(1:K,ER_FOOBI,...
 'marker','^','markersize',makerSize,...
   'linestyle','none','color',black_o,'LineWidth',LineWidth);
d42 = semilogy(1:1,ER_FOOBI(1),...
    'marker','^','markersize',makerSize,...
    'linestyle','-','color',black_o,'LineWidth',LineWidth);

d0 = semilogy(1:K,ER_TenSOFO,...
    'linestyle','-','color',red_o,'LineWidth',LineWidth);
d01 = plot(1:K,ER_TenSOFO,...
 'marker','s','markersize',makerSize,...
   'linestyle','none','color',red_o,'LineWidth',LineWidth);
d02 = semilogy(1:1,ER_TenSOFO(1),...
    'marker','s','markersize',makerSize,...
    'linestyle','-','color',red_o,'LineWidth',LineWidth);

lgd = legend([ d22,  d42,  d02], '\texttt{SOBIUM}',...
      '\texttt{FOOBI}', '\texttt{TenSOFO}');
lgd.FontSize = 18;
set(lgd, 'Interpreter', 'latex', 'Color', [0.95, 0.95, 0.95]);
xlabel('SNR (dB)','interpreter','latex','FontSize',13,'FontName','Times New Roman');
ylabel('$\texttt{RE}(\mathbf{A},\mathbf{A}_{\mathrm{est}})$','interpreter','latex','FontSize',13,'FontName','Times New Roman');

% 
h=gca;
set(gca,'YScale', 'log');
set(h,'FontSize',16,'XGrid','on','YGrid','on','GridLineStyle','-','MinorGridLineStyle',':','FontName','Times New Roman');
xticks([1 2 3 4 5])
xticklabels({'0','5','10','15','20'})
yticks([1e-2 1e-1 1])
yticklabels({'10^{-2}','10^{-1}','10^{0}'})
set(h,'FontSize', 30);
axis([1 5.5 0.5e-2 1e0])
grid on;
box on;





