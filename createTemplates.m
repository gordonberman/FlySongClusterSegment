function [outputData,allPeakIdx,allNormalizedPeaks,peakAmplitudes,isNoise,scores,options] = ...
                                createTemplates(data,options,plotsOn)
    
    %Inputs:
                 
    %data -> 1-D array of data points
    %options -> self-explanatory
    %plotsOn -> plot final templates? [T (default)/ F]

    
    %Outputs:
    
    %templates -> L x 1 cell array containing M_i x d arrays of template peaks
    %allPeakIdx -> all found peak locations
    %allNormalizedPeaks -> N x d array containing all found peaks 
    %                       corresponding to the locations in allPeakIdx
    %peakAmplitudes -> N x 1 vector of peak normalization factors
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
    
    %finds the noise level for the data set 
    if options.noiseLevel < 0
        
        fprintf(1,'   Finding Noise Level\n');
        if length(data) > 1e6
            shortdata = data(1:1e6);
        else
            shortdata = data;
        end
        
        out = returnShuffledPhaseData(shortdata);
        options.noiseLevel = std(out);
        
    end

    
    fprintf(1,'   Finding Preliminary Peak Locations\n');
    [normalizedPeaks,peakIdx,~,peakAmplitudes] = ...
        findNormalizedPeaks(data,options.noiseLevel,options.sigmaThreshold,...
        options.diffThreshold,options.smoothSigma,[],options.minNoiseLevel);
    
    
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
    amplitudes = cell(options.k,1);
    for i=1:options.k
        templates{i} = normalizedPeaks(idx==i,:);
        amplitudes{i} = peakAmplitudes(idx==i);
    end
    
    
          
    splitted = true;
    isNoise = false(size(templates));
    while splitted
        
        noiseTemplates = templates(isNoise);
        noiseAmplitudes = amplitudes(isNoise);
        isNoiseOld = isNoise(isNoise);
        
        [templates2,isNoise,splitted,amplitudes2] = selectTemplates(templates(~isNoise),amplitudes,false);

        templates = [templates2;noiseTemplates];
        amplitudes = [amplitudes2;noiseAmplitudes];
        
        isNoise = [isNoise;isNoiseOld];
        
    end
    
    
    templateSizes = zeros(length(templates),1);
    for i=1:length(templates)
        templateSizes(i) = length(templates{i}(:,1));
    end
    templates = templates(templateSizes >= options.min_template_size);
    isNoise = isNoise(templateSizes >= options.min_template_size);
    amplitudes = amplitudes(templateSizes >= options.min_template_size);
    
    outputData.isNoise = isNoise;
    outputData.templates = templates;
    outputData.amplitudes = amplitudes;   
    
    close all
    
    if plotsOn
        makeTemplateHistograms(templates,histogramBins);
    end
    
    
    L = length(templates);
    d = length(templates{1}(1,:));
    coeffs = cell(L,1);
    projStds = cell(L,1);
    means = cell(L,1);
    L_templates = zeros(L,1);
    fprintf(1,'   Finding Template Bases and Projections\n');
    for i=1:L
        
        fprintf(1,'      Template #%2i\n',i);
        L_templates(i) = length(templates{i}(:,1));
        
        %find Data Set Mean
        means{i} = mean(templates{i});
    
        %perform PCA on set of normalized peaks
        [coeffs{i},scores,~] = princomp(templates{i});
        
        projStds{i} = std(scores);
    end
    
    
    %adjust bases sets for sub-sampled data sets
    minLength = min(L_templates);
    if minLength < 2*d
        q = round(minLength / 2);
    else
        q = d;
    end
    
    for i=1:L
        coeffs{i} = coeffs{i}(:,options.first_mode:q);
        projStds{i} = projStds{i}(options.first_mode:q);
    end
    
    %find baseline noise levels
    fprintf(1,'   Calculating Baseline Noise Levels\n');
    if options.use_likelihood_threshold
        baselines = repmat(options.baseline_threshold,sum(~isNoise),1);
    else
        [baselines,~] = findTemplateBaselines(templates,coeffs,means,projStds,isNoise,options.baseline_quantile);
    end
    
    outputData.coeffs = coeffs;
    outputData.projStds = projStds;
    outputData.L_templates = L_templates;
    outputData.means = means;
    outputData.baselines = baselines;

    isNoise = outputData.isNoise;
    
    
    