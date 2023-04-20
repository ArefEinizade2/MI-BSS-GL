function [S, Lambda] = GL_From_Covvectors_vech_RealData(V, Lambda_vec)
N = size(V,1);
%%
M_p = full(duplication_matrix(N));
vech_idx = pinv(M_p)*vec(eye(N));
Sel = eye(length(vech_idx));
Sel = Sel(vech_idx==0,:);
M_N = pinv(Sel);
M = M_p*M_N;
%%
for i = 1:N
   V_p(:,i) = vec(V(:,i)*V(:,i)'); 
end
A = pinv(M)*V_p;
A_pinv = pinv(A);
%%
m = N*(N-1)/2;
s = ones(m,1);
w = ones(m,1);
Lambda = Lambda_vec(:);
z1 = s - A*Lambda;
z2 = w - s;
k = 0;
rho1 = 1;
rho2 = 1;
cnt1 = 1000;
cnt2 = 1000;
e = 1e-2;
%%
s_k = s; s_k_plus = s_k;
w_k = w; w_k_plus = w_k;
L_k = Lambda; L_k_plus = L_k;
z1_k = z1; z1_k_plus = z1_k;
z2_k = z2; z2_k_plus = z2_k;
%%
Start = true;
MaxIter = 10;
t_k_plus = 1;
t_k = 1;
% while ( (Start) ||  ( (norm(s_k_plus-s_k)/norm(s_k))>e &&  (norm(w_k_plus-w_k)/norm(w_k))>e && (norm(L_k_plus-L_k)/norm(L_k))>e && ...
%         (norm(z1_k_plus-z1_k)/norm(z1_k))>e && (norm(z2_k_plus-z2_k)/norm(z2_k))>e ) )
% x = (rho1*A*L_k + rho2*w_k + z1_k - z2_k)/(rho1 + rho2);
% for iter = 1:MaxIter
while ( (Start) ||  norm(L_k_plus-L_k)/norm(L_k)>e )
    %% initialization:
    Start = false;
    k = k + 1;
    s_k = s_k_plus;
    w_k = w_k_plus;
    L_k = L_k_plus;
    z1_k = z1_k_plus;
    z2_k = z2_k_plus;
%     t_k = t_k_plus;
    %% Update S:
%     param.maxit = 10000;
%     param.nu = 1e-8;
    param.verbose = 0;

    x = (rho1*A*L_k + rho2*w_k + z1_k - z2_k)/(rho1 + rho2);
    gamma = (1/(rho1 + rho2));
    s_k_plus = prox_l1(x, gamma, param);

%     t_k_plus = (1+sqrt(1+4*t_k^2))/2;
%     x = s_k_plus + ((t_k-1)/t_k_plus)*(s_k_plus - s_k);
    
    %% Update Lambda:
    L_k_plus = A_pinv*(s_k_plus - z1_k/rho1);
    %% Update W:
    w_k_plus = s_k_plus + z2_k/rho2;
    w_k_plus(w_k_plus<0) = 0;
%     for i = 1:100
%         P_w_k_plus = w_k_plus;
%         P_w_k_plus(w_k_plus<0) = 0;
%         w_k_plus = w_k_plus + 0.7*(P_w_k_plus - w_k_plus);
%         
%     end
    %% Update z1 & z2:
    z1_k_plus = z1_k - (s_k_plus - A*L_k_plus)*rho1;
    z2_k_plus = z2_k - (w_k_plus - s_k_plus)*rho2;
    %%
    rho1 = rho1 * cnt1;
    rho2 = rho2 * cnt2;
end
%%
s = abs(A*L_k_plus);
S = reshape(M*s,[N,N]); S = S - diag(diag(S));
% if norm(S,'fro')~=0
%     S = S / max(S(:));
% end
% [~,d] = eig(S);
Lambda = L_k_plus;
end