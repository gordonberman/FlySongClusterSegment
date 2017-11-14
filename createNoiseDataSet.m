function [noiseData,fs,psd,maxP,sigmaNoise,threshold,percentNoise] = ...
    createNoiseDataSet(data,noiseDataLength,options)
%Creates a noise "data set" from the array 'data'
%
%Inputs:
%   data -> array of data or a string pointing to a .wav file 
%   noiseDataLength -> length of returned noise data
%   options -> options structure (see makeParameterStructure.m for details)
%
%Outputs:
%   noiseData -> array of noise data
%   fs -> frequencies in noise data power spectrum
%   psd -> power spectral density from noise data
%   maxP -> maximum psd value (note: not the frequency)
%   sigmaNoise -> standard deviation of noise
%   threshold -> threshold for determining signal from noise
%   percentNoise -> fraction of 'data' classified as noise


    addpath(genpath('./utilities/'));
    addpath(genpath('./subroutines/'));
    
    
    if nargin < 3 || isempty(options)
        options.setAll = true;
    else
        options.setAll = false;
    end
    options = makeParameterStructure(options);
    
    
    if ischar(data)
        [data,Fs] = audioread(data);
        options.Fs = Fs;
    end
    
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
    
    
    y = log10(gaussianfilterdata(data.^2,sigma));
    obj = findBestGMM_AIC(y,maxNumGaussians,replicates,maxNumPoints);
    
    idx = argmax(obj.mu);
    minIdx = argmin(obj.mu);
    posts = posterior(obj,y);
    posts = posts(:,idx);
    
    mask = (posts >= noise_posterior_threshold | y > obj.mu(idx)) & y > obj.mu(minIdx);
    threshold = min(y(mask));
    
    if threshold < min_noise_threshold
        mask = y > min_noise_threshold;
    end
    
    CC = bwconncomp(~mask);
    numSegments = CC.NumObjects;
    
    d = [(1:numSegments)' returnCellLengths(CC.PixelIdxList)'];
    idx = randperm(numSegments);
    
    d = d(idx,:);
    totalLengths = sum(d(:,2));
    if totalLengths < noiseDataLength
        noiseDataLength = totalLengths;
        idx2 = numSegments;
    else
        idx2 = find(cumsum(d(:,2)) >= noiseDataLength,1,'first');
    end

    
    noiseData = zeros(noiseDataLength,1);
    count = 0;
    for i=1:idx2
        noiseData(count + (1:d(i,2))) = data(CC.PixelIdxList{d(i,1)});
        count = count + d(i,2);
    end
        
        
    fs = linspace(0,Fs/2,ceil(length(noiseData)/2));
    f = fft(noiseData);
    p = f.*conj(f);
    p = p(2:length(fs));
    psd = p ./ (sum(p)*(fs(2)-fs(1)));
    fs = fs(2:end);
    
    maxP = max(psd);

  
    sigmaNoise = std(noiseData);
    percentNoise = mean(mask);
        
    
        
        
        
    
