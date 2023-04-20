function Beta = SFD(y)
% Implements the Score_Function_Difference (SFD) for two estimated sources as proposed in the following paper:
% Babaie-Zadeh, Massoud, and Christian Jutten. "A general approach for mutual information minimization and its application to blind source separation." Signal Processing 85, no. 5 (2005): 975-995.
%% In the case of finding its contents useful for your research work, kindly please also cite our paper addressed below:
% [1] Einizade, Aref, and Sepideh Hajipour Sardouie. "Joint Graph Learning and Blind Separation of Smooth Graph Signals Using Minimization of Mutual Information and Laplacian Quadratic Forms." IEEE Transactions on Signal and Information Processing over Networks 9 (2023): 35-47.
%%
% Usage: 
%   >> Beta = SFD(y);
%   y: nxT, n is the number of sources and must be 2, T is the number of temporal samples
%   Beta: nxT: the Score_Function_Difference

%% 
Psi = MSF(y);

Phi = JSF(y);

Beta = Psi - Phi;

end