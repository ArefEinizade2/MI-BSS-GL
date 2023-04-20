function SNR = getSNR_BSS(s, y)
% A function to calculate the SNR output metric in source separation performance

S = normalize(s ,2,'range');

P = size(s, 1);

Order = 1 : P;

combinations = perms(Order);

SNR_comb = zeros(1, size(combinations, 1));

for i = 1 : size(combinations,1)
    
    y_i = normalize(y(combinations(i, :), :) ,2,'range');
    SNR_comb(i) = mean(10 * log10 (mean(S.^2, 2) ./  mean( ( S - y_i ).^2, 2))); 
    
end

SNR = max(SNR_comb);

end