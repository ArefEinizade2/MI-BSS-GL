function [ZjadeGL, B] = GraphJADEGL_KG(X, P, N, DelayNum, Lambda_b, Maxiter, W_est)

% known graphs:

WinNum = floor(size(X,2)/N);
X = X(:,1:WinNum*N);
Graph_decorrelation_term = Lambda_b;
%% ICA_Preprocessing:
[WightedX, Whitening_mat] = ICA_Preprocessing(X);
%% C matrices:
[C_cell] = Compute_Comulants(WightedX);
%%

    %% S matrices:
S_cell = zeros(P, (P*P)*DelayNum);
%     S_cell = zeros(P, WinNum*(P*P)*DelayNum);
for i = 1:WinNum
    X_tild_i = WightedX(:,(i-1)*N+1:i*N);
    S_cell_WinNum = zeros(P, DelayNum*P*P);
    for k = 1:DelayNum
        SS_mat = zeros(P,P*P);
        for p = 1:P
            SS_mat(:,(p-1)*P+1:p*P) = S_matrice_RealData(X_tild_i, W_est{p}, k);
        end
            S_cell_WinNum(:, (k-1)*(P*P)+1:(k)*(P*P)) = SS_mat;
    end
%         S_cell(:, (i-1)*((P*P)*DelayNum)+1:i*((P*P)*DelayNum)) = S_cell_WinNum;
    S_cell = S_cell + S_cell_WinNum;
end
S_cell = S_cell/WinNum;
%%
%     if iter ~= 1
%         value1 = Compute_GraphJADE_critria(S_cell);
%         value2 = Compute_GraphJADE_critria(C_cell);
%         Max = value2/(value1+value2);
%         Graph_decorrelation_term = Max;%Graph_decorrelation_term + (Max/(1*Maxiter));
%     end

%% Joint Diagonalziation:
B = Joint_diagonalization_jader_version([(Graph_decorrelation_term)*S_cell, (1-Graph_decorrelation_term)*C_cell], WightedX, Whitening_mat);

ZjadeGL = B * X;


end

