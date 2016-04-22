function [templates,allPeakIdx,allNormalizedPeaks,isNoise,scores,options] = createTemplates(data,options,plotsOn)

    
    %Inputs:
                 
    %data -> 1-D array of data points
    %options -> self-explanatory
    %plotsOn -> plot final templates? [T (default)/ F]

    
    %Outputs:
    
    %templates -> L x 1 cell array containing M_i x d arrays of template peaks
    %allPeakIdx -> all found peak locations
    %allNormalizedPeaks -> N x d array containing all found peaks 
    %                       corresponding to the locations in allPeakIdx
    %isNoise -> L x 1 binary array, returns zero if marked signal, 1 if noise
    %scores -> peak pca projections
    %options -> returned options structure

    addpath(genpath('./chronux'))

    histogramBins = 100;
    
        
    if nargin < 2 || isempty(options)
        options.setAll = true;
    else
        options.setAll = false;
    end
    options = makeDefaultOptions(options);

    if nargin < 3 || isempty(plotsOn)
        plotsOn = true;
    end
    
    
%     if options.noiseLevel <= 0
%         fprintf(1,'   Finding Noise Level\n');
%         if numel(data) > 1e6
%             shortdata = data(1:1e6);
%         else
%             shortdata = data;
%         end
%         [ssf] = sinesongfinder(shortdata,options.fs,12,20,.1,.01,.05,[0 1000]); %returns ssf, which is structure containing the following fields
%         noise = findnoise(ssf,3,80,1000);
%         options.noiseLevel = noise.sigma;
%     end

    options.noiseLevel = quantile(data,.95);
    options.sigmaThreshold = 1;

    
    fprintf(1,'   Finding Preliminary Peak Locations\n');
    [normalizedPeaks,peakIdx,~] = ...
        findNormalizedPeaks(data,options.noiseLevel,options.sigmaThreshold,...
        options.diffThreshold,options.smoothSigma);
    
    
    [~,scores,~] = princomp(normalizedPeaks);
    
    N = length(peakIdx);
    if N > options.maxNumPeaks
        q = randperm(N);
        q = q(1:options.maxNumPeaks);
        normalizedPeaks = normalizedPeaks(q,:);
        peakIdx = peakIdx(q);
        scores = scores(q,:);
    end
        
    allNormalizedPeaks = normalizedPeaks;
    allPeakIdx = peakIdx;
    scores = scores(:,1:options.template_pca_dimension);
    
    
    fprintf(1,'   Clustering Peaks\n');
    kmeans_options = statset('MaxIter',options.kmeans_maxIter);
    idx = kmeans(scores,options.k,'replicates',options.kmeans_replicates,'options',kmeans_options);
        
    templates = cell(options.k,1);
    for i=1:options.k
        templates{i} = normalizedPeaks(idx==i,:);
    end
    
    
          
    splitted = true;
    isNoise = false(size(templates));
    while splitted
        
        noiseTemplates = templates(isNoise);
        isNoiseOld = isNoise(isNoise);
        
        [templates2,isNoise,splitted] = selectTemplates(templates(~isNoise),false);

        templates = [templates2;noiseTemplates];
        
        isNoise = [isNoise;isNoiseOld];
        
    end
    
    
    templateSizes = zeros(length(templates),1);
    for i=1:length(templates)
        templateSizes(i) = length(templates{i}(:,1));
    end
    templates = templates(templateSizes >= options.min_template_size);
    isNoise = isNoise(templateSizes >= options.min_template_size);
    
    
    close all
    
    if plotsOn
        makeTemplateHistograms(templates,histogramBins);
    end
    
    