function S = UtilFunc3(X_tild, W, k)

    P = size(X_tild,1);
    
    S = P*(X_tild*(W^k)*X_tild')/norm((W^k)*X_tild','fro');

end
