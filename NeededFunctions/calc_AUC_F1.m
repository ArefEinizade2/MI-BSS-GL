function [AUC, F1] = calc_AUC_F1(W_org, W_est)

% A function to calculate the AUC and F1 metrics in graph recovery performance

[~,~,~,AUC] = perfcurve(vec(full(squareform_sp(W_org))), vec(full(squareform_sp(W_est))),1);

[tpr, fpr] = roc(vec(full(squareform_sp(W_org)))', vec(full(squareform_sp(W_est)))');

F = (tpr) ./ (fpr + tpr); F = F(~isnan(F)); F1 = mean(F);

end