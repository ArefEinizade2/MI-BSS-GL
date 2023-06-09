function [y, B] = MI_BSS_KG(x, L, param)
% Implements the MI-BSS-KG method for blind separation of smooth graph signal sources with known graphs, as proposed in our paper.
% In the case of finding its contents useful for your research work, kindly please also cite our paper addressed below:
% [1] Einizade, Aref, and Sepideh Hajipour Sardouie. "Joint Graph Learning and Blind Separation of Smooth Graph Signals Using Minimization of Mutual Information and Laplacian Quadratic Forms." IEEE Transactions on Signal and Information Processing over Networks 9 (2023): 35-47.
%%
% Usage: 
%   >> [y, B] = MI_BSS_KG(x, param);

%   x: nxT, n is the number of sources, T is the number of observed temporal samples
%   param.mu: the learning rate used in the Gradient Descent step, e.g., 0.1
%   param.lambda: balances the weight between minimization of mutual information and graph smoothness terms, e.g., 1e-2
%   param.WinNum: the number of windows in the cyclostationary sources, e.g., 40
%   L: a n-length cell containg the original Laplacian matrices of the graph signal sources, where L{i} has the size of T/WinNum x T/WinNum
%   param.mu: the learning rate used in the Gradient Descent step, e.g., 0.1
%   param.lambda: balances the weight between minimization of mutual information and graph learning terms, e.g., 1e-2
%   param.WinNum: the number of windows in the cyclostationary sources, e.g., 40
%   param.Tol: tolerance for reaching convergence
%   param.MaxIter: maximum iterations for convergence

%   y: nxT, the estimated sources
%   B: nxn, the estimated unmixing matrix
%% Initialization:
mu = param.mu;
lambda = param.lambda;
WinNum = param.WinNum;
Tol = param.Tol;
MaxIter = param.MaxIter;

[P, M] = size(x);

N = M / WinNum;

B = eye(size(x,1));

%%  Loop for blind separation of smooth graph signal sources with Known graphs:  

for iter = 1 : 5000

    disp(['MIBSS-GL-KG iter ', num2str(iter)])

    y = B * x;

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

    MI_grad_reg = zeros(size(MI_grad));

    for p = 1 : P

        Phi_p = zeros(size(B));

        for r = 1 : WinNum

            Phi_p = Phi_p + x(:, (r - 1) * N + 1 : r * N) * L{p} * x(:, (r - 1) * N + 1 : r * N)';

        end

        Phi_p = Phi_p / WinNum;

        MI_grad_reg(p, :) = B(p,:) * Phi_p';

    end

    MI_grad = MI_grad + ((lambda) * (MI_grad_reg));

    % Gradient Descent Step:
    new_B = B - mu * MI_grad;

    sigma = std(y, [], 2);

    new_B = new_B ./ sigma;

    error(iter) = norm((new_B - B), 'fro')/norm(B,'fro');      

    B = new_B;

    if error(iter) < Tol           
        break;
    end

end

end