function [WightedX, B] = ICA_Preprocessing(X)
%   B = jadeR(X, m) is an m*n matrix such that Y=B*X are separated sources
%    extracted from the n*T data matrix X.
%   If m is omitted,  B=jadeR(X)  is a square n*n matrix (as many sources as sensors)
%
% Blind separation of real signals with JADE.  Version 1.9.   August 2013
%
% Usage: 
%   * If X is an nxT data matrix (n sensors, T samples) then
%     B=jadeR(X) is a nxn separating matrix such that S=B*X is an nxT
%     matrix of estimated source signals.
%   * If B=jadeR(X,m), then B has size mxn so that only m sources are
%     extracted.  This is done by restricting the operation of jadeR
%     to the m first principal components. 
%   * Also, the rows of B are ordered such that the columns of pinv(B)
%     are in order of decreasing norm; this has the effect that the
%     `most energetically significant' components appear first in the
%     rows of S=B*X.
%
% Quick notes (more at the end of this file)
%
%  o this code is for REAL-valued signals.  An implementation of JADE
%    for both real and complex signals is also available from
%    http://perso.telecom-paristech.fr/~cardoso/guidesepsou.html
%
%  o This algorithm differs from the first released implementations of
%    JADE in that it has been optimized to deal more efficiently
%    1) with real signals (as opposed to complex)
%    2) with the case when the ICA model does not necessarily hold.
%
%  o There is a practical limit to the number of independent
%    components that can be extracted with this implementation.  Note
%    that the first step of JADE amounts to a PCA with dimensionality
%    reduction from n to m (which defaults to n).  In practice m
%    cannot be `very large' (more than 40, 50, 60... depending on
%    available memory and CPU time)
%
%  o See more notes, references and revision history at the end of
%    this file and more stuff on the WEB
%    http://perso.telecom-paristech.fr/~cardoso/guidesepsou.html
%
%  o This code is supposed to do a good job!  Please report any
%    problem to cardoso@iap.fr
[n,T]	= size(X);
m = n;
X	= X - mean(X')' * ones(1,T);
[U,D]     = eig((X*X')/T) ; %% An eigen basis for the sample covariance matrix
[Ds,k]    = sort(diag(D)) ; %% Sort by increasing variance
PCs       = n:-1:n-m+1    ; %% The m most significant princip. comp. by decreasing variance

%% --- PCA  ----------------------------------------------------------
B         = U(:,k(PCs))'    ; % At this stage, B does the PCA on m components

%% --- Scaling  ------------------------------------------------------
scales    = sqrt(Ds(PCs)) ; % The scales of the principal components .
B         = diag(1./scales)*B  ; % Now, B does PCA followed by a rescaling = sphering


%% --- Sphering ------------------------------------------------------
X         = B*X;  %% We have done the easy part: B is a whitening matrix and X is white.
WightedX = X;
clear U D Ds k PCs scales ;

end
