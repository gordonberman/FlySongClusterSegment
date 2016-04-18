function [baselines,projectionLikes] = findTemplateBaselines(templates,coeffs,means,projStds,isNoise,quantileValue)


    L_signal = sum(~isNoise);
    L_noise = sum(isNoise);
    baselines = zeros(L_signal,1);
    noiseVals = find(isNoise);
    
    if quantileValue < 0
        
        baselines(:) = -9e99;
        projectionLikes = [];
        
    else
        
        projectionLikes = cell(L_signal,1);
        for i=1:L_signal
            projectionLikes{i} = zeros(length(templates{i}(:,1)),L_noise);
            for k=1:L_noise
                j = noiseVals(k);
                projections = (templates{i} - repmat(means{j},length(templates{i}(:,1)),1))*coeffs{j};
                projectionLikes{i}(:,k) = findDataSetLikelihoods(projections,zeros(size(means{i})),projStds{j});
            end
            
            baselines(i) = quantile(projectionLikes{i}(:),quantileValue);
            
        end
        
    end