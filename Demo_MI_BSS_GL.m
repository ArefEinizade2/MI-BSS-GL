%% This a demo for simultaneous graph learning and blind separation of smooth graph signal sources
%% In the case of finding its contents useful for your research work, kindly please cite the following paper. Thanks a lot for your attention.
% [1] Einizade, Aref, and Sepideh Hajipour Sardouie. "Joint Graph Learning and Blind Separation of Smooth Graph Signals Using Minimization of Mutual Information and Laplacian Quadratic Forms." IEEE Transactions on Signal and Information Processing over Networks 9 (2023): 35-47.
clc; clear; close all;    
%% It takes about 120 minutues to run for 100 independent realizations and P={2,3,4} sources
%% Generates the results of Figures 1 and 2, i.e., the subsection named "BSS and GL performance" in our manuscript
%% Add necessary functions:
rng(2);
addpath('./GraphJADE_UtilFunctions'); 
addpath('./GraphJADEGL_UtilFunctions'); 
addpath('./FDPG'); % the exploited graph learning approach from the following paper:
addpath('./NeededFunctions');
%% setting input parameters:
P_rd = 0.6; % edge probability of Erdos-Renyi graphs
N = 10; % number of nodes of Erdos-Renyi graphs
param_ER.connected = 1; % make the generated Erdos-Renyi graphs to be connected 
MeanVec = zeros(1, N); % mean of the generated graph signals
DelayNum = 1; % Delay parameter in the GraphJADE and GraphJADE-GL methods
param.mu = 0.1; % the learning rate in the gradient descent step of the MI-BSS-GL and MI-BSS-KG methods
param.lambda = 1e-2; % lambda in the MI-BSS-GL method
param.WinNum = 40; % number of windows, omega in the paper, for the MI-BSS-GL method
param.Tol = 1e-3; % tolerance for reaching convergence
param.MaxIter = 5000; %  maximum iterations for convergence
nSamples = 1; % number of i.i.d. graph signals per window
b_GraphJADE = 0.4; % balancing parameter in the GraphJADE and GraphJADE-GL methods
Maxiter_GraphJADE = 10; % number of convergence iterations of the GraphJADE and GraphJADE-GL methods
MaxGen = 100; % number of genertated independent realizations 
sigma_vec = [0.3, 0.5, 0.7, 0.9]; % the span of variation of the noise level (sigma) 
source_num_vec = 2:4;  % the span of variation of the source numbers 

%% Perform the related analysis explained in our manuscript:
t0 = tic;
for source_num = source_num_vec

    AUC = [];
    
    F1 = [];
    
for sigma = sigma_vec
    
%%
for iter = 1 : MaxGen
    
    source_num
    sigma
    iter

s = zeros(source_num, param.WinNum * N);

G = {}; W_org = {}; L = {}; V = {}; d = {}; D = {}; gftcoeff = {};

for source = 1 : source_num
    
    G{source} = gsp_erdos_renyi(N, P_rd, param_ER); % generating Erdos-Renyi graphs using GSPBOX

    W_org{source} = full(G{source}.W); % original adjacency matrix

    L{source} = diag(sum(W_org{source}, 2)) - W_org{source}; % original laplacian matrix

    [V{source}, D{source}] = eig(L{source}); % EVD of laplacian matrix

    d{source} = pinv(D{source}); % pseudo-inverse of eigenvalue matrix of the laplacian

    for r = 1 : param.WinNum % generate smooth graph signals for each window
        
        gftcoeff{source} = mvnrnd(MeanVec, d{source}, nSamples);

        s(source, (r - 1) * N + 1 : r * N ) = V{source} * gftcoeff{source}' + sigma * randn(size(V{source} * gftcoeff{source}'));
        
    end
end

s = mapstd(s);

%% Mix sourses and generate observations
A = 0.8 * rand(size(s,1), size(s,1));
A = A - diag(diag(A));
A = A + eye(size(s,1));
X = A * s;
%% Apply MI-BSS method (MI-BSS-GL with lambda=0): 
param_MIBSS = param;
param_MIBSS.lambda = 0;
[y, B1] = MI_BSS_GL(X, param_MIBSS);

% calculate the Minimum Distance criterion metric in source sepration:
MD1(find(source_num_vec==source_num), find(sigma_vec==sigma), iter) = N * (source_num - 1) * Minimum_Distance_crit(A, B1);

% calculate the SNR output metric in source separation performance:
S = normalize(s ,2,'range');
y = normalize(y ,2,'range');
SNR_Final_MIBSS(find(source_num_vec==source_num), find(sigma_vec==sigma), iter) = mean(10 * log10 (mean(S.^2,2) ./  mean((S-y).^2 , 2)));  % Equivalent
%% Apply the proposed MI-BSS-GL method:
[y_GS, B2, W_est] = MI_BSS_GL(X, param);

% calculate the Minimum Distance criterion metric in source sepration:
MD2(find(source_num_vec==source_num), find(sigma_vec==sigma), iter) = N * (source_num - 1) * Minimum_Distance_crit(A, B2);

% calculate the SNR output metric in source separation performance:
S = normalize(s ,2,'range');
y_GS = normalize(y_GS ,2,'range');
SNR_Final_MIBSS_GL(find(source_num_vec==source_num), find(sigma_vec==sigma), iter) = mean(10 * log10 (mean(S.^2,2) ./  mean((S-y_GS).^2 , 2)));  % Equivalent

% calculate the AUC and F1 metrics in graph recovery performance:
for p = 1 : size(s,1)
    
    [AUC(find(sigma_vec==sigma), p, iter), F1(find(sigma_vec==sigma), p, iter)] = calc_AUC_F1(W_org{p}, W_est{p});

end
%% Apply MI-BSS-KG method as gold standard:

[y_GS_KG, B2_KG] = MI_BSS_KG(X, L, param);

% calculate the Minimum Distance criterion metric in source sepration:
MD2_KG(find(source_num_vec==source_num), find(sigma_vec==sigma), iter) = N * (source_num - 1) * Minimum_Distance_crit(A, B2_KG);

% calculate the SNR output metric in source separation performance:
SNR_Final_MIBSS_GS_KG(find(source_num_vec==source_num), find(sigma_vec==sigma), iter) = getSNR_BSS(s, y_GS_KG);

%% Apply the JADE method:
B3 =  jadeRR(X, source_num);

y3 = B3 * X;

% calculate the Minimum Distance criterion metric in source sepration:
MD3(find(source_num_vec==source_num), find(sigma_vec==sigma), iter) = N * (source_num - 1) * Minimum_Distance_crit(A, B3);

% calculate the SNR output metric in source separation performance:
SNR_Final_JADE(find(source_num_vec==source_num), find(sigma_vec==sigma), iter) = getSNR_BSS(s, y3);
%% Apply the GraphJADE-GL method:

[~, B4] = GraphJADEGL(X, source_num, N, DelayNum, b_GraphJADE, Maxiter_GraphJADE);

% calculate the Minimum Distance criterion metric in source sepration:
MD4(find(source_num_vec==source_num), find(sigma_vec==sigma), iter) = N * (source_num - 1) * Minimum_Distance_crit(A, B4);

y4 = B4 * X;

% calculate the SNR output metric in source separation performance:
SNR_Final_GraphJADE(find(source_num_vec==source_num), find(sigma_vec==sigma), iter) = getSNR_BSS(s, y4);

end


end

AUC_cell{find(source_num_vec==source_num)} = AUC;

F1_cell{find(source_num_vec==source_num)} = F1;


end
%% Plot MD error results:

figure;
for source_num = source_num_vec
    
    for sigma = sigma_vec

        MD1_vec(find(sigma_vec==sigma)) = round(mean(squeeze(MD1(find(source_num_vec==source_num), find(sigma_vec==sigma), :))),3);
        
        MD2_vec(find(sigma_vec==sigma)) = round(mean(squeeze(MD2(find(source_num_vec==source_num), find(sigma_vec==sigma), :))),3);

        MD2_KG_vec(find(sigma_vec==sigma)) = round(mean(squeeze(MD2_KG(find(source_num_vec==source_num), find(sigma_vec==sigma), :))),3);
        
        MD3_vec(find(sigma_vec==sigma)) = round(mean(squeeze(MD3(find(source_num_vec==source_num), find(sigma_vec==sigma), :))),3);

        MD4_vec(find(sigma_vec==sigma)) = round(mean(squeeze(MD4(find(source_num_vec==source_num), find(sigma_vec==sigma), :))),3);

        xx = squeeze(MD2(find(source_num_vec==source_num), find(sigma_vec==sigma), :)) - squeeze(MD1(find(source_num_vec==source_num),...
            find(sigma_vec==sigma), :)); %% Statistical test
        
        [h1,MDMI_p1(find(source_num_vec==source_num), find(sigma_vec==sigma)),ci1,stats1] = ttest(xx);

    end
        subplot(1, length(source_num_vec), find(source_num_vec==source_num)); 

        plot(sigma_vec, MD1_vec, '-*', 'LineWidth', 2); hold on; 
        
        xlabel('\sigma'); ylabel('N(P-1)ave(D^2)'); title(['P: ', num2str(source_num)])

        plot(sigma_vec, MD2_vec, '-O', 'LineWidth',2);
        
        plot(sigma_vec, MD2_KG_vec, '-O', 'LineWidth',2);
        
        plot(sigma_vec, MD3_vec, '-*', 'LineWidth',2);

        plot(sigma_vec, MD4_vec, '-x', 'LineWidth',2);
        
        xticks(sigma_vec); xlim([min(sigma_vec(:)), max(sigma_vec(:))])
        

end

legend(['MI-BSS'], ['MI-BSS-GL'], ['MI-BSS-KG'], ['JADE'], ['GraphJADEGL'],...
    'Location', 'best', 'Orientation','horizontal')

%% Plot the output SNRs:
figure;
for source_num = source_num_vec
    
    for sigma = sigma_vec

        SNR_Final_1_vec(find(sigma_vec==sigma)) = mean(SNR_Final_MIBSS(find(source_num_vec==source_num), find(sigma_vec==sigma), :), 'all');
        
        SNR_Final_2_vec(find(sigma_vec==sigma)) = mean(SNR_Final_MIBSS_GL(find(source_num_vec==source_num), find(sigma_vec==sigma), :), 'all');
        
        SNR_Final_3_vec(find(sigma_vec==sigma)) = mean(SNR_Final_MIBSS_GS_KG(find(source_num_vec==source_num), find(sigma_vec==sigma), :), 'all');

        SNR_Final_4_vec(find(sigma_vec==sigma)) = mean(SNR_Final_JADE(find(source_num_vec==source_num), find(sigma_vec==sigma), :), 'all');

        SNR_Final_5_vec(find(sigma_vec==sigma)) = mean(SNR_Final_GraphJADE(find(source_num_vec==source_num), find(sigma_vec==sigma), :), 'all');
    
    end
    
        subplot(1, length(source_num_vec), find(source_num_vec==source_num)); 

        plot(sigma_vec, SNR_Final_1_vec, '-*', 'LineWidth', 2); hold on; xlabel('\sigma', 'FontWeight', 'bold', 'FontSize', 14); 
        
        ylabel('SNR (dB)', 'FontWeight', 'bold', 'FontSize', 14); title(['P: ', num2str(source_num)], 'FontWeight', 'bold', 'FontSize', 16)

        plot(sigma_vec, SNR_Final_2_vec, '-O', 'LineWidth',2);
        
        plot(sigma_vec, SNR_Final_3_vec, '-*', 'LineWidth',2);

        plot(sigma_vec, SNR_Final_4_vec, '-O', 'LineWidth',2);

        plot(sigma_vec, SNR_Final_5_vec, '-*', 'LineWidth',2);

        xticks(sigma_vec); xlim([min(sigma_vec(:)), max(sigma_vec(:))])
        
        legend(['MI-BSS'], ['MI-BSS-GL'], ['MI-BSS-KG'], ['JADE'], ['GraphJADE-GL'], 'Location', 'best')

end
%% Plot the AUC and F1-measure metrics for graph learning performance:

figure;

subplot(1, 2, 1); 

for source_num = source_num_vec
    
        AAA = AUC_cell{find(source_num_vec==source_num)};
    
        AAA = mean(AAA, [2,3]);

        plot(sigma_vec, AAA, '-O', 'LineWidth',2); hold on;
        
        title('AUC vs. Noise level (\sigma)'); ylabel('AUC'); xlabel('\sigma');

        xticks(sigma_vec); xlim([min(sigma_vec(:)), max(sigma_vec(:))])

end

legend('P=2', 'P=3', 'P=4', 'Location', 'bestoutside', 'Orientation','horizontal');

% F1:

subplot(1, 2, 2); 

for source_num = source_num_vec
    
        AAA = F1_cell{find(source_num_vec==source_num)};
    
        AAA = mean(AAA, [2,3]);
        
        plot(sigma_vec, AAA, '-O', 'LineWidth',2); hold on;
        
        title('F1 vs. Noise level (\sigma)'); ylabel('F1'); xlabel('\sigma');

        xticks(sigma_vec); xlim([min(sigma_vec(:)), max(sigma_vec(:))])

end
legend('P=2', 'P=3', 'P=4', 'Location', 'bestoutside', 'Orientation','horizontal');

%%
t1 = toc(t0);

disp(['>>>>> run-time: ', num2str(round(t1/60,2)), ' minutes']);

