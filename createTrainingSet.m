function [trainingSet,thresholds] = createTrainingSet(data,trainingSetLength,minFromEachDataSet,options)


    addpath(genpath('./utilities/'));
    addpath(genpath('./subroutines/'));
    
    if ~iscell(data)
        data = {data};
    end
    
    N = length(data);
    
    masks = cell(N,1);
    thresholds = zeros(N,1);
    numSegments = zeros(N,1);
    
    if nargin < 3 || isempty(minFromEachDataSet)
        minFromEachDataSet = trainingSetLength / (N*10);
    end
    
    if nargin < 4 || isempty(options)
        options.setAll = true;
    else
        options.setAll = false;
    end
    options = makeParameterStructure(options);
    
    
    Fs = options.fs;
    sigma = 10*options.smoothingLength_noise * Fs / 1000;
    maxNumGaussians = options.maxNumGaussians_noise;
    replicates = options.replicates_GMM;
    maxNumPoints = options.maxNumPeaks;
    noise_posterior_threshold = options.noise_posterior_threshold;
    if options.min_noise_threshold > 0
        min_noise_threshold = log10(options.min_noise_threshold.^2);
    else
        min_noise_threshold = [];
    end
    high_pass_filter_cutoff = options.high_pass_filter_cutoff / (Fs/2);
    butterworth_order = options.butterworth_order;
    if high_pass_filter_cutoff > 0
        [b,a] = butter(butterworth_order,high_pass_filter_cutoff,'high');
    end
    
    
    CCs = cell(N,1);
    for i=1:N
        
        if high_pass_filter_cutoff > 0
            data{i} = filter(b,a,data{i});
        end
        
        y = log10(gaussianfilterdata(data{i}.^2,sigma));
        obj = findBestGMM_AIC(y,maxNumGaussians,replicates,maxNumPoints);
        
        idx = argmax(obj.mu);
        minIdx = argmin(obj.mu);
        posts = posterior(obj,y);
        posts = posts(:,idx);
        
        masks{i} = (posts >= noise_posterior_threshold | y > obj.mu(idx)) & y > obj.mu(minIdx);
        thresholds(i) = min(y(masks{i}));
        
        if thresholds(i) < min_noise_threshold
            thresholds(i) = min_noise_threshold;
            masks{i} = y > min_noise_threshold;
        end
       
        CCs{i} = bwconncomp(masks{i});
        numSegments(i) = CCs{i}.NumObjects;
        
    end
        
    d = zeros(sum(numSegments),3);
    count = 0;
    for i=1:N
        d((1:numSegments(i)) + count,1) = i;
        d((1:numSegments(i)) + count,2) = 1:numSegments(i);
        d((1:numSegments(i)) + count,3) = returnCellLengths(CCs{i}.PixelIdxList);
        count = count + numSegments(i);
    end
        
    dataSetMins = zeros(N,1);
    for i=1:N
        dataSetMins(i) = min([sum(d(d(:,1)==i,3)),minFromEachDataSet]);
    end
    
    
    idx = randperm(sum(numSegments));
    d = d(idx,:);
    
    cumSums = cumsum(d(:,3));
    idx = find(cumSums >= trainingSetLength,1,'first');
    if isempty(idx)
        idx = sum(numSegments);
%     else
%         numFromDataSets = zeros(N,1);
%         for i=1:N
%             numFromDataSets(i) = sum(d(d(1:idx,1)==i,3));
%         end
%         idx2 = find(numFromDataSets < dataSetMins);
%         if ~isempty(idx2)
%             
%             for i=1:length(idx2)
%                 
%                 
%                 
%             end
%             
%         end
        
    end
    
    
    
    trainingSet = zeros(sum(d(1:idx,3)),1);
    trainingSet_origins = zeros(idx,3);
    count = 0;
    for i=1:idx
        trainingSet((1:d(i,3)) + count) = data{d(i,1)}(CCs{d(i,1)}.PixelIdxList{d(i,2)});
        trainingSet_origins(i,:) = [count+1,count+d(i,3),d(i,1)];  
        count = count + d(i,3);
    end
    
    
    
    
    
    
    
    
    
    