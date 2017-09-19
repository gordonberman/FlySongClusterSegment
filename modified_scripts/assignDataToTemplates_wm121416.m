function [groupings,peakIdxGroup,likes,allPeakIdx,allNormalizedPeaks,noiseThreshold] = ...
                     assignDataToTemplates(data,outputData,options)

                 
    %Inputs:
                 
    %data -> 1-D array of data points
    %outputData -> structure containing information about templates
    %options -> self-explanatory

    
    %Outputs:
    
    %groupings -> L x 1 cell array, each cell containing an N_i x d array 
    %               of grouped peaks
    %peakIdxGroup -> L x 1 cell array of index values for each 
    %               corresponding group
    %likes -> N x L array of the log-likelihood values that each peak 
    %           belongs to the corresponding model
    %allPeakIdx -> all found peak locations
    %allNormalizedPeaks -> N x d array containing all found peaks 
    %                       corresponding to the locations in allPeakIdx
    %coeffs -> L x 1 cell array of model PCA bases
    %projStds -> L x 1 cell array of model projection standard deviations
    
    
    addpath(genpath('./utilities/'));
    addpath(genpath('./subroutines/'));
    
    templates = outputData.templates;
    projStds = outputData.projStds;
    coeffs = outputData.coeffs;
    means = outputData.means;
    
    L = length(templates);
    d = length(templates{1}(1,:));
       
    if nargin < 3 || isempty(options)
        options.setAll = true;
    else
        options.setAll = false;
    end
    options = makeParameterStructure(options);
    
    Fs = options.fs;

    
    %find all potential peaks in the new data set
    fprintf(1,'   Finding Peak Locations\n');

    diffThreshold = d;

    maxNumGaussians_noise = options.maxNumGaussians_noise;
    maxNumPeaks_GMM = options.maxNumPeaks;
    replicates_GMM = options.replicates_GMM; 
    smoothingLength_noise = options.smoothingLength_noise * Fs / 1000;
    minRegionLength = round(options.minRegionLength * Fs / 1000);
    min_noise_threshold = options.min_noise_threshold;
    median_filter_length = round(options.median_filter_length * Fs / 1000);
    noise_posterior_threshold = options.noise_posterior_threshold;
    diff_threshold_multiplier = options.diff_threshold_multiplier;
    if min_noise_threshold > 0
        min_noise_threshold = log10(min_noise_threshold.^2);
    else
        min_noise_threshold = [];
    end
    
    high_pass_filter_cutoff = options.high_pass_filter_cutoff / (Fs/2);
    butterworth_order = options.butterworth_order;
    if high_pass_filter_cutoff > 0
        [b,a] = butter(butterworth_order,high_pass_filter_cutoff,'high');
        data = filter(b,a,data);
    end
    
    
    if median_filter_length > 0
        if mod(median_filter_length,2) == 0
            median_filter_length = median_filter_length + 1;
        end
        data = medfilt1(data,median_filter_length);
    end
    
    
    [newData,~,noiseThreshold,peakIdx] = filterDataAmplitudes_wm121416(data,...
        smoothingLength_noise,minRegionLength,maxNumGaussians_noise,...
        replicates_GMM,maxNumPeaks_GMM,diffThreshold*diff_threshold_multiplier,[],...
        min_noise_threshold,noise_posterior_threshold);
    
    
    r = (diffThreshold-1)/2;
    peakIdx = peakIdx(peakIdx > r & peakIdx < length(newData)-r); %added to eliminate edge cases
    N = length(peakIdx);
    normalizedPeaks = zeros(N,diffThreshold);
    peakAmplitudes = zeros(N,1);
    for i=1:N
        a = newData(peakIdx(i) + (-r:r));
        peakAmplitudes(i) = sqrt(mean(a.^2));
        normalizedPeaks(i,:) = a./peakAmplitudes(i).*sign(newData(peakIdx(i)));
    end
    
    allPeakIdx = peakIdx;
    allNormalizedPeaks = normalizedPeaks;
    
            
    %find likelihood values
    fprintf(1,'   Finding Likelihoods\n');
    likes = zeros(N,L);
    
    for i=1:L
        
        fprintf(1,'      Template #%2i\n',i);
        
        projections = bsxfun(@minus,normalizedPeaks,means{i})*coeffs{i};
                
        likes(:,i) = findDataSetLikelihoods(projections,zeros(size(means{i})),projStds{i});
    end
    
    clear projections 
    
    
    %assign peaks to the template with maximum likelihood
    [~,idx] = max(likes,[],2);
    
    
    
    groupings = cell(L,1);
    peakIdxGroup = cell(L,1);
    for i=1:L
        q = idx == i;
        groupings{i} = allNormalizedPeaks(q,:);
        peakIdxGroup{i} = allPeakIdx(q);
    end
    
    if isfield(outputData,'templateGroupings') && ~isempty(outputData.templateGroupings)
        
        k = unique(outputData.templateGroupings);
        M = length(k);
        groupings2 = cell(M,1);
        peakIdxGroup2 = cell(M,1);

        for i=1:M
            groupings2{i} = cell2mat(groupings(outputData.templateGroupings==k(i)));
            peakIdxGroup2{i} = cell2mat(peakIdxGroup(outputData.templateGroupings==k(i)));
        end
        
        groupings = groupings2;
        peakIdxGroup = peakIdxGroup2;
        
    end
    
    
    
    
    