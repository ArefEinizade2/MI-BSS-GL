function S = S_matrice_RealData(X_tild, W, k)

    N = size(W,1);
    P = size(X_tild,1);
    
    if norm(W-zeros(size(W)),'fro')~=0       
        S = P*(X_tild*(W^k)*X_tild')/norm(W^k*X_tild','fro');
%         S = P*(X_tild*(W^k)*X_tild')/norm(X_tild*W^k*X_tild','fro'); dge b parameterha hasas nemishe injoori!
        
    else
        S = ((X_tild)*(W^k)*(X_tild'))/(N-k);
    end

end