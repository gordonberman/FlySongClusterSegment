function outputData = clusterTemplates_wm112116(outputData,numClusters,options,plotsOn)

    if nargin < 3 
        options = [];
    end
    options = makeParameterStructure(options);
    
    if nargin < 4 || isempty(plotsOn)
        plotsOn = true;
    end
        
    isNoise = outputData.isNoise;
    
    
    [~,~,~,~,~,D] = ...
        calculateTemplateFrequencyProfiles(outputData.templates,options.fs);

    z = makePdistForm(D);
    Y = linkage(z,'single');
    a = cluster(Y,'maxclust',numClusters);
    
    k = unique(a);
    if a(end) ~= max(k)
        b = a;
        b(a==a(end)) = length(k);
        idx = setdiff(k,a(end));
        for i=1:length(idx)
            b(a == k(idx(i))) = i;
        end
        a = b;
    end
    outputData.templateGroupings = a;
    
    k = unique(a);
    outputData.isNoiseTemplateGrouping = false(length(k),1);
    for i=1:length(k)
        outputData.isNoiseTemplateGrouping(i) = min(isNoise(a==k(i)));
    end

    
    if plotsOn
        
        histogramBins = 100;
        k = unique(a);
        for i=1:length(k)
            figure
            makeTemplateHistograms_wm112116(outputData.templates(a==k(i)),histogramBins,...
                [],[-.5 .5]*sqrt(outputData.diffThreshold),...
                outputData.isNoise(a==k(i)),outputData.percentBelowNoiseThreshold(a==k(i))); %modified 11/17/16
        end
        
    end
    
    
    
    
    