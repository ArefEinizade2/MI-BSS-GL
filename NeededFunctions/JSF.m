function Phi = JSF(y)
% Implements the Joint Score Function (JSF) for two estimated sources as proposed in the following paper:
% Babaie-Zadeh, Massoud, and Christian Jutten. "A general approach for mutual information minimization and its application to blind source separation." Signal Processing 85, no. 5 (2005): 975-995.
%% In the case of finding its contents useful for your research work, kindly please also cite our paper addressed below:
% [1] Einizade, Aref, and Sepideh Hajipour Sardouie. "Joint Graph Learning and Blind Separation of Smooth Graph Signals Using Minimization of Mutual Information and Laplacian Quadratic Forms." IEEE Transactions on Signal and Information Processing over Networks 9 (2023): 35-47.
%%
% Usage: 
%   >> Phi = JSF(y);
%   y: nxT, n is the number of sources and must be 2, T is the number of temporal samples
%   Phi: nxT: the Joint Score Function
%%
Moment = @(k,i) mean(y(i,:).^k); 

Moment_Cross = @(k1, k2, i1, i2) mean((y(i1,:).^k1) .* (y(i2,:).^k2));

Phi = zeros(2, size(y,2)); 

E = [1 Moment(1,1), Moment(2,1), Moment(3,1), Moment(1,2), Moment(2,2), Moment(3,2);
    Moment(1,1), Moment(2,1), Moment(3,1), Moment(4,1), Moment_Cross(1,1,1,2), Moment_Cross(1,2,1,2), Moment_Cross(1,3,1,2); 
    Moment(2,1), Moment(3,1), Moment(4,1), Moment(5,1), Moment_Cross(2,1,1,2), Moment_Cross(2,2,1,2), Moment_Cross(2,3,1,2); 
    Moment(3,1), Moment(4,1), Moment(5,1), Moment(6,1), Moment_Cross(3,1,1,2), Moment_Cross(3,2,1,2), Moment_Cross(3,3,1,2);
    Moment(1,2), Moment_Cross(1,1,2,1), Moment_Cross(1,2,2,1), Moment_Cross(1,3,2,1), Moment(2,2), Moment(3,2), Moment(4,2);
    Moment(2,2), Moment_Cross(2,1,2,1), Moment_Cross(2,2,2,1), Moment_Cross(2,3,2,1), Moment(3,2), Moment(4,2), Moment(5,2);
    Moment(3,2), Moment_Cross(3,1,2,1), Moment_Cross(3,2,2,1), Moment_Cross(3,3,2,1), Moment(4,2), Moment(5,2), Moment(6,2)];

M1 = [0; 1; 2*Moment(1,1); 3*Moment(2,1); 0; 0; 0];

M2 = [0; 0; 0; 0; 1; 2*Moment(1,2); 3*Moment(2,2)];
            
theta1 = E \ M1;

theta2 = E \ M2;

Phi(1,:) = theta1(1) + theta1(2)*y(1,:) + theta1(3)*y(1,:).^2 + theta1(4)*y(1,:).^3 + ...
    theta1(5)*y(2,:) + theta1(6)*y(2,:).^2 + theta1(7)*y(2,:).^3;
    
Phi(2,:) = theta2(1) + theta2(2)*y(1,:) + theta2(3)*y(1,:).^2 + theta2(4)*y(1,:).^3 + ...
    theta2(5)*y(2,:) + theta2(6)*y(2,:).^2 + theta2(7)*y(2,:).^3;

end