function [y, B, W_est] = MI_BSS_GL(x, mu, lambda, WinNum)
% Implements the MI-BSS-GL method for simultaneous graph learning and blind separation of smooth graph signal sources, as proposed in our paper.
% In the case of finding its contents useful for your research work, kindly please also cite our paper addressed below:
% [1] Einizade, Aref, and Sepideh Hajipour Sardouie. "Joint Graph Learning and Blind Separation of Smooth Graph Signals Using Minimization of Mutual Information and Laplacian Quadratic Forms." IEEE Transactions on Signal and Information Processing over Networks 9 (2023): 35-47.
%%
% Usage: 
%   >> [y, B, W_est] = MI_BSS_GL(x, mu, lambda, WinNum);

%   x: nxT, n is the number of sources, T is the number of observed temporal samples
%   mu: the learning rate used in the Gradient Descent step, e.g., 0.1
%   lambda: balances the weight between minimization of mutual information and graph learning terms, e.g., 1e-2
%   WinNum: the number of windows in the cyclostationary sources, e.g., 40
%   s: nxT, the originial independent sources for claculating the quality of sepration performance

%   y: nxT, the estimated sources
%   B: nxn, the estimated unmixing matrix
%   W_est: a n-length cell containg the learned adjacency matrices of the graph signal sources, where W_est{i} has the size of T/WinNum x T/WinNum

%% Initialization:

    [P, M] = size(x);
    
    N = M / WinNum;
    
    B = eye(size(x,1));

    MI_grad_reg = zeros(size(B)); error_g = 0;
              
%%  Loop for simultaneous graph learning and blind separation of smooth graph signal sources:  
    
    for iter = 1 : 5000
        
        disp(['MIBSS-GL iter ', num2str(iter)])
                        
        y = B * x;
        
        %% >>>>>>>>>>>> Graph Learning >>>>>>>>>>>>>>>>>>>>>>>
        if lambda ~= 0
            for p = 1 : P

                Y_p = reshape(y(p, :), [N, WinNum]);

                Z_p = sparse(gsp_distanz(Y_p').^2);

                z_p = squareform_sp(Z_p/WinNum);

                W_est{p} = FDPG_log_degree(z_p, 1e0, 1e0, 1e-3); % learned adjacency matrix

                L{p} = diag(sum(W_est{p}, 2)) - W_est{p}; % learned Laplacian matrix
                
            end

        %% >>>>>>>>>>>> Blind Source Separartion using the minimization of grapg regularized mutual information >>>>>>>>>>>>>>>>>>>>>>>
            MI_grad_reg = zeros(size(B));

            for p = 1 : P

                Phi_p = zeros(size(B));

                for r = 1 : WinNum

                    Phi_p = Phi_p + x(:, (r - 1) * N + 1 : r * N) * L{p} * x(:, (r - 1) * N + 1 : r * N)';

                end

                Phi_p = Phi_p / WinNum;

                MI_grad_reg(p, :) = B(p,:) * Phi_p';

            end
        end
        
        Gamma = zeros(size(x));
        % Pairwise estimation of Score Difference Functions:
        for i = 1:P-1
                        
            for j = i+1:P
                                                
                Beta_star = SFD([y(i,:); y(j,:)]);

                Beta = Beta_star;
                
                Gamma(i,:) = Gamma(i,:) + Beta(1,:);
                
                Gamma(j,:) = Gamma(j,:) + Beta(2,:);

            end
                        
        end
                
        MI_grad = Gamma*x'/size(x,2);
        
        % Gradient Descent Step:
        
        MI_grad = MI_grad + ((lambda) * (MI_grad_reg));
        
        new_B = B - mu * MI_grad;
                        
        sigma = std(y, [], 2);
                
        new_B = new_B ./ sigma;
                
        error(iter) = norm((new_B - B), 'fro')/norm(B,'fro');      
        
        B = new_B;
                                        
        if error(iter) < 1e-3            
            break;
        end
        
    end
    
    error_g = mean(error_g(:, end), 'all');
            
end