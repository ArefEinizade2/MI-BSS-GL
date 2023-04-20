function [MD_square, MD] = Minimum_Distance_crit(ActualMixingMatrix, EstimatedUnmixingMatrix)
% A function to calculate the Minimum Distance criterion metric in source
% separation performance, as proposed in the following reference:
% [1] Ilmonen, P., Nordhausen, K., Oja, H., & Ollila, E. (2010). A new performance index for ICA: properties, computation and asymptotic analysis. In Latent Variable Analysis and Signal Separation: 9th International Conference, LVA/ICA 2010, St. Malo, France, September 27-30, 2010. Proceedings 9 (pp. 229-236). Springer Berlin Heidelberg.

    P = size(ActualMixingMatrix,1);
    Per = perms(1:P);
    I = eye(P);
    K = size(Per,1);
    criterion = zeros(1,K);
    G = EstimatedUnmixingMatrix*ActualMixingMatrix;
    G_New = G;
    for i = 1:size(G,1)
        for j = 1:size(G,2)          
            G_New(i,j) = (G(i,j)^2)/(norm(G(i,:))^2 + 1e-20);           
        end
    end
    
    for i = 1:K        
        C = I(:,Per(i,:));       
        criterion(i) = trace(C*G_New);
    end
    
    MD_square = abs((P-max(criterion)))/(P-1);
    MD = sqrt(MD_square);

end