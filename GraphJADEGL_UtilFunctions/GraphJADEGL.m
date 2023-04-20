function [ZjadeGL, U, W_est] = GraphJADEGL(X, P, N, DelayNum, Lambda_b, Maxiter)
% Implements the GraphJADEGL method for simultaneous graph learning and blind separation of graph signal sources, as proposed in our papers [1,2].
% In the case of finding its contents useful for your research work, kindly please also cite our paper addressed below:
% [1] Einizade, Aref, Sepideh Hajipour Sardouie, and Mohammad B. Shamsollahi. "Simultaneous graph learning and blind separation of graph signal sources." IEEE Signal Processing Letters 28 (2021): 1495-1499.
% [2]Einizade, Aref, and Sepideh Hajipour Sardouie. "A unified approach for simultaneous graph learning and blind separation of graph signal sources." IEEE Transactions on Signal and Information Processing over Networks 8 (2022): 543-555.
%%
WinNum = floor(size(X,2)/N);
X = X(:,1:WinNum*N);
Graph_decorrelation_term = 0;
U = eye(P);
%% Optimization Process:
for p = 1:P
W_init = zeros(N,N); W_init = (W_init+W_init')/2; W_init = W_init - diag(diag(W_init));  
W_est{p} = W_init;
end
%% ICA_Preprocessing:
[WightedX, Whitening_mat] = ICA_Preprocessing(X);
%% C matrices:
[C_cell] = Compute_Comulants(WightedX);
%%
for iter = 1:Maxiter
    disp(['GraphJADE-GL iter #',num2str(iter)])

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
    U = Joint_diagonalization_jader_version([(Graph_decorrelation_term)*S_cell, (1-Graph_decorrelation_term)*C_cell], WightedX, Whitening_mat);
    ZjadeGL = U*X;
    Graph_decorrelation_term = Lambda_b;
    %% Covariance Estimation and Optimize subject to Adjacency matrices and their eigenvalues:
        for p = 1:P
            Z_p = ZjadeGL(p,:);
            Z_p_mat = reshape(Z_p,[length(Z_p)/WinNum,WinNum]);
            Cov_est = cov(Z_p_mat');
            [V,D] = eig(Cov_est);
            [W_est{p}, GraphEigenValues{p}] = GL_From_Covvectors_vech_RealData(V, diag(D));
            ComponentsEigenValues{p} = diag(D);
            W = W_est{p};
%             W_est{p} = W/max(W(:));
                      
        end
        
%         if norm(U - U_new, 'fro')/norm(U, 'fro') < 1e-3
%             break;
%         end
%         U = U_new;

end


end

