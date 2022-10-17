function [Axon,Stray] = merge_clusters(DF,IDs,Posterior,Weighting,Round)

Stray = false(size(DF,1),1);
Axon = zeros(length(unique(IDs)),size(DF,2));

for Ax = 1:length(unique(IDs))
    TempDistance = sort(Posterior(IDs == Ax,:)');
    TempRatio = TempDistance(end,:) - TempDistance(end-1,:);
    
    % don't cluster when one posterior is at least 95% of that at
    % second closest cluster
    if any(TempRatio<=0.95) && Round == 1
        Stray(find(IDs == Ax)) = true;
    else
        
        TempBoutons = DF(IDs == Ax,:);
        
        % weighted average by skewness
        if strcmp(Weighting,'Skewness')
            TempSkew = skewness(TempBoutons,[],2);
            TempSkew(TempSkew>5) = 5;
            TempSkew = (TempSkew - mean(TempSkew)) + 1;
            for C = 1:length(TempSkew); TempBoutons(C,:) = TempBoutons(C,:) .* TempSkew(C);end
        end
        
        Axon(Ax,:) = nanmean(TempBoutons,1);
    end
end

