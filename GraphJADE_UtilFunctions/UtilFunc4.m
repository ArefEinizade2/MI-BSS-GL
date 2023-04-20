function U = UtilFunc4(CM, WightedX, B )
%%% Joint diagonalization of the cumulant matrices used in JADE method:
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T = size(WightedX,2);
m = size(CM,1);
nbcm = size(CM,2)/m;
%% Init
if 0, 	%% Init by diagonalizing a *single* cumulant matrix.  It seems to save
	% some computation time `sometimes'.  Not clear if initialization is really worth
	% it since Jacobi rotations are very efficient.  On the other hand, it does not
	% cost much...

	[V,D]	= eig(CM(:,1:m)); % Selectng a particular cumulant matrix.
	for u=1:m:m*nbcm,         % Accordingly updating the cumulant set given the init
		CM(:,u:u+m-1) = CM(:,u:u+m-1)*V ; 
	end;
	CM	= V'*CM;

else,	%% The dont-try-to-be-smart init
	V	= eye(m) ; % la rotation initiale
end;

%% Computing the initial value of the contrast 
Diag    = zeros(m,1) ;
On      = 0 ;
Range   = 1:m ;
for im = 1:nbcm,
  Diag  = diag(CM(:,Range)) ;
  On    = On + sum(Diag.*Diag) ;
  Range = Range + m ;
end
Off = sum(sum(CM.*CM)) - On ;


seuil	= 1.0e-6 / sqrt(T) ; % A statistically scaled threshold on `small' angles
encore	= 1;
sweep	= 0; % sweep number
updates = 0; % Total number of rotations
upds    = 0; % Number of rotations in a given seep
g	= zeros(2,nbcm);
gg	= zeros(2,2);
G	= zeros(2,2);
c	= 0 ;
s 	= 0 ;
ton	= 0 ;
toff	= 0 ;
theta	= 0 ;
Gain    = 0 ;

%% Joint diagonalization proper

while encore, encore=0;   

  sweep = sweep+1;
  upds  = 0 ; 
  Vkeep = V ;
  
  for p=1:m-1,
    for q=p+1:m,

      Ip = p:m:m*nbcm ;
      Iq = q:m:m*nbcm ;
      
      %%% computation of Givens angle
      g	    = [ CM(p,Ip)-CM(q,Iq) ; CM(p,Iq)+CM(q,Ip) ];
      gg    = g*g';
      ton   = gg(1,1)-gg(2,2); 
      toff  = gg(1,2)+gg(2,1);
      theta = 0.5*atan2( toff , ton+sqrt(ton*ton+toff*toff) );
      Gain  = (sqrt(ton*ton+toff*toff) - ton) / 4 ;
      
      %% Givens update
      if abs(theta) > seuil,
%      if Gain > 1.0e-3*On/m/m ,
	encore  = 1 ;
	upds    = upds    + 1;
	c	= cos(theta); 
	s	= sin(theta);
	G	= [ c -s ; s c ] ;
	
	pair 		= [p;q] ;
	V(:,pair) 	= V(:,pair)*G ;
	CM(pair,:)	= G' * CM(pair,:) ;
	CM(:,[Ip Iq]) 	= [ c*CM(:,Ip)+s*CM(:,Iq) -s*CM(:,Ip)+c*CM(:,Iq) ] ;
	

	On   = On  + Gain;
	Off  = Off - Gain;
	
	%% fprintf('jade -> %3d %3d %12.8f\n',p,q,Off/On);
      end%%of the if
    end%%of the loop on q
  end%%of the loop on p
  updates = updates + upds ;
  
end%%of the while loop


%%% A separating matrix
%   ===================
B	= V'*B ;


%%% Permut the rows of the separating matrix B to get the most energetic components first.
%%% Here the **signals** are normalized to unit variance.  Therefore, the sort is
%%% according to the norm of the columns of A = pinv(B)

A           = pinv(B) ;
[Ds,keys]   = sort(sum(A.*A)) ;
B           = B(keys,:)       ;
B           = B(m:-1:1,:)     ; % Is this smart ?


% Signs are fixed by forcing the first column of B to have non-negative entries.

b	= B(:,1) ;
signs	= sign(sign(b)+0.1) ; % just a trick to deal with sign=0
B	= diag(signs)*B ;
U = B;
end