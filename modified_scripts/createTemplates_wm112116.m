function [outputData,allPeakIdx,allNormalizedPeaks,peakAmplitudes,isNoise,allScores,options] = ...
                                createTemplates_wm112116(data,options,plotsOn)
    
    %Inputs:
                 
    %data -> 1-D array of data points
    %options -> self-explanatory (default = [])
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

    
    addpath(genpath('./utilities/'));
    addpath(genpath('./subroutines/'));

    histogramBins = 100;
        
    if nargin < 2 || isempty(options)
        options.setAll = true;
    else
        options.setAll = false;
    end
    options = makeParameterStructure(options);

    Fs = options.fs;  
    
    
    if nargin < 3 || isempty(plotsOn)
        plotsOn = true;
    end
    
    fprintf(1,'   Finding Peak Locations\n');
    
    %load run parameters
    maxNumGaussians_noise = options.maxNumGaussians_noise;
    maxNumPeaks_GMM = options.maxNumPeaks;
    maxNumPeaks = options.maxNumPeaks;
    replicates_GMM = options.replicates_GMM; 
    smoothingLength_noise = options.smoothingLength_noise * Fs / 1000;
    minRegionLength = round(options.minRegionLength * Fs / 1000);
    maxIPI = options.maxIPI / 1000;
    IPI_sigma = options.IPI_sigma / 1000; %* Fs / 1000; 
    num_IPI_halfWidths = options.num_IPI_halfWidths;
    amplitude_threshold = options.amplitude_threshold; 
    numIPIBins = 10000;
    min_noise_threshold = options.min_noise_threshold;
    median_filter_length = round(options.median_filter_length * Fs / 1000);
    min_seperation = floor(options.min_seperation * Fs / 1000);
    noise_posterior_threshold = options.noise_posterior_threshold;
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
    
    
    
    
    %find regions masked by pulse amplitude
    [~,~,noiseThreshold,smoothedPeakLocations] = ...
        filterDataAmplitudes(data,smoothingLength_noise,minRegionLength,...
                        maxNumGaussians_noise,replicates_GMM,...
                        maxNumPeaks_GMM,min_seperation,[],...
                        min_noise_threshold,noise_posterior_threshold);
     
    %find IPIs                
    IPIs = diff(smoothedPeakLocations)/options.fs;
    IPIs = IPIs(IPIs < maxIPI);
    
    %find IPI distribution through kernel density estimation
    [Y,X] = hist(IPIs,numIPIBins);
    Y = normalizeHist(X,Y);
    IPI_sigma = IPI_sigma / (X(2) - X(1));
    Y = gaussianfilterdata(Y,IPI_sigma);
    s = fit(X',Y','spline');
    
    
    %estimate mode of the IPI distribution
    modeIPI_estimate = fminsearch(@(x) -s(x),X(argmax(Y)));
    peakValue = s(modeIPI_estimate);
    
    %estimate width of IPI distribution
    halfMax_location = fzero(@(x) s(x)-peakValue/2,X(argmax(Y)));
    halfMax_location = modeIPI_estimate - ...
        num_IPI_halfWidths*abs(modeIPI_estimate-halfMax_location);
    min_location = fminsearch(@(x) s(x),halfMax_location);
    threshold_location = max(halfMax_location,min_location);
    diffThreshold = floor(threshold_location*Fs);
    diffThreshold = max([diffThreshold 3*smoothingLength_noise]);
    
    %make sure that diffThreshold is an odd number of points
    if mod(diffThreshold,2) == 0
        diffThreshold = diffThreshold - 1;
    end
    
    %refilter data set using the found noise threshold and diffThreshold
    [newData,~,~,peakIdx] = filterDataAmplitudes(data,...
        smoothingLength_noise,minRegionLength,[],[],[],...
        diffThreshold,noiseThreshold,[],noise_posterior_threshold);
    
    %find normalized peaks
    %N = length(peakIdx);
    r = (diffThreshold-1)/2;
    peakIdx = peakIdx(peakIdx > r & peakIdx < length(newData)-r); %added to eliminate edge cases
    N = length(peakIdx);
    normalizedPeaks = zeros(N,diffThreshold);
    peakAmplitudes = zeros(N,1);
    signs = sign(newData(peakIdx));
    
    for i=1:N
        a = newData(peakIdx(i) + (-r:r));
        peakAmplitudes(i) = sqrt(mean(a.^2));
        normalizedPeaks(i,:) = (a./peakAmplitudes(i)).*signs(i);
    end

    %save output data
    outputData.noiseThreshold = noiseThreshold;
    outputData.diffThreshold = diffThreshold;
    outputData.min_location = min_location;
    outputData.halfMax_location = halfMax_location;
    outputData.modeIPI_estimate = modeIPI_estimate;
    outputData.IPI_width = abs(modeIPI_estimate-halfMax_location);
    if options.template_pca_dimension > min(diffThreshold,N-1)
        options.template_pca_dimension = min(diffThreshold,N-1);
    end
    
    [~,scores,~] = pca(normalizedPeaks);
    
    
    if N > maxNumPeaks
        q = randperm(N,options.maxNumPeaks);
        normalizedPeaks = normalizedPeaks(q,:);
        peakIdx = peakIdx(q);
        scores = scores(q,:);
    end
        

    allNormalizedPeaks = normalizedPeaks;
    allPeakIdx = peakIdx;
    allScores = scores(:,1:options.template_pca_dimension);
    scores = allScores;
    
    fprintf(1,'   Clustering Peaks\n');
    kmeans_options = statset('MaxIter',options.kmeans_maxIter);
           
    idx = kmeans(scores,options.k,'replicates',options.kmeans_replicates,'options',kmeans_options);
    clusterIdx = idx;
        
    templates = cell(options.k,1);
    amplitudes = cell(options.k,1);
    for i=1:options.k
        templates{i} = normalizedPeaks(idx==i,:);
        amplitudes{i} = peakAmplitudes(idx==i);
    end
    templates = templates(returnCellLengths(amplitudes) > 0);
    amplitudes = amplitudes(returnCellLengths(amplitudes) > 0);
    
    isNoise = false(size(templates));
    
    
    %uncomment this to allow for human annotation of templates
    %     splitted = true;
    %     isNoise = false(size(templates));
    %     while splitted
    %
    %         noiseTemplates = templates(isNoise);
    %         noiseAmplitudes = amplitudes(isNoise);
    %         isNoiseOld = isNoise(isNoise);
    %
    %         [templates2,isNoise,splitted,amplitudes2] = selectTemplates(templates(~isNoise),amplitudes,false);
    %
    %         templates = [templates2;noiseTemplates];
    %         amplitudes = [amplitudes2;noiseAmplitudes];
    %
    %         isNoise = [isNoise;isNoiseOld];
    %
    %     end
    
    
    %automatically determine noise templates
    percentBelowNoiseThreshold = zeros(size(templates));
    noiseVal = log10(sqrt(10.^noiseThreshold));
    for i=1:length(templates)
        percentBelowNoiseThreshold(i) = mean(log10(amplitudes{i}) < noiseVal);
    end 
    isNoise = percentBelowNoiseThreshold > amplitude_threshold | isNoise;
    
    
    [~,idx] = sort(double(isNoise));
    templates = templates(idx);
    amplitudes = amplitudes(idx);
    isNoise = isNoise(idx);
    percentBelowNoiseThreshold = percentBelowNoiseThreshold(idx);
    outputData.percentBelowNoiseThreshold = percentBelowNoiseThreshold;
        

    %First, error out if no templates are generated
    if isempty(templates)
        error('     No templates generated!');
    end
    
    templateSizes = zeros(length(templates),1);
    for i=1:length(templates)
        templateSizes(i) = length(templates{i}(:,1));
    end
    templates = templates(templateSizes >= options.min_template_size);
    isNoise = isNoise(templateSizes >= options.min_template_size);
    amplitudes = amplitudes(templateSizes >= options.min_template_size);
    percentBelowNoiseThreshold = percentBelowNoiseThreshold(templateSizes >= options.min_template_size); %modified 11/17/16
    
    outputData.isNoise = isNoise;
    outputData.templates = templates;
    outputData.amplitudes = amplitudes;  
    outputData.percentBelowNoiseThreshold = percentBelowNoiseThreshold; %modified 11/17/16
    
    close all
    
    %make template plots
    if plotsOn
        figure
        makeTemplateHistograms_wm112116(templates,histogramBins,[],...
            [-.5 .5]*sqrt(outputData.diffThreshold),isNoise,...
            percentBelowNoiseThreshold); %modified 11/17/16
    end
    
    
    %defining templates
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
        [coeffs{i},scores,~] = pca(templates{i});
        
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
    
    
    outputData.coeffs = coeffs;
    outputData.projStds = projStds;
    outputData.L_templates = L_templates;
    outputData.means = means;

    isNoise = outputData.isNoise;
    
    
    if options.run_tsne
        fprintf(1,'   Computing t-SNE Embedding\n');
        options.signalLabels = [];
        options.tsne_readout = 25;
        [yData,~,~,~] = run_tSne(allScores,options);
        outputData.yData = yData;
        outputData.idx = idx;
        figure
        scatter(yData(:,1),yData(:,2),[],clusterIdx,'filled')
        colormap(jet)
        axis equal tight off 
        colorbar
        set(gca,'fontsize',16,'fontweight','bold')
    end
    
    