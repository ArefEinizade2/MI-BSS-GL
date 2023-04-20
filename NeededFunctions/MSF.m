function Psi = MSF(y)
% Implements the Marginal Score Function (MSF) for two estimated sources as proposed in the following paper:
% Babaie-Zadeh, Massoud, and Christian Jutten. "A general approach for mutual information minimization and its application to blind source separation." Signal Processing 85, no. 5 (2005): 975-995.
%% In the case of finding its contents useful for your research work, kindly please also cite our paper addressed below:
% [1] Einizade, Aref, and Sepideh Hajipour Sardouie. "Joint Graph Learning and Blind Separation of Smooth Graph Signals Using Minimization of Mutual Information and Laplacian Quadratic Forms." IEEE Transactions on Signal and Information Processing over Networks 9 (2023): 35-47.
%%
% Usage: 
%   >> Psi = jader(y);
%   y: nxT, n is the number of sources and must be 2, T is the number of temporal samples
%   Psi: nxT: the Marginal Score Function
%%

Moment = @(k,i) mean(y(i,:).^k);

Psi = zeros(2, size(y,2));

for i = 1 : 2 
    
    E = [1 Moment(1,i), Moment(2,i), Moment(3,i);
        Moment(1,i), Moment(2,i), Moment(3,i), Moment(4,i);
        Moment(2,i), Moment(3,i), Moment(4,i), Moment(5,i);
        Moment(3,i), Moment(4,i), Moment(5,i), Moment(6,i)];
    
    M = [0; 1; 2*Moment(1,i); 3*Moment(2,i)];
    
    theta = E \ M;

    Psi(i,:) = theta(1) + theta(2)*y(i,:) + theta(3)*y(i,:).^2 + theta(4)*y(i,:).^3;
    
end

end