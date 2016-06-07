function [normalizedPeaks,peakIdx,data,amplitudes] = ...
    findNormalizedPeaks(data,s,threshold,diffThreshold,smoothSigma,maxNumber,minNoiseLevel)

    
    if nargin < 3 || isempty(threshold)
        threshold = 5;
    end
    
    if nargin < 4 || isempty(diffThreshold)
        diffThreshold = 150;
    end
    
    if nargin < 5 || isempty(smoothSigma)
        smoothSigma = 5;
    end
    
    if nargin < 7 || isempty(minNoiseLevel)
        minNoiseLevel = 0;
    end
    
    if smoothSigma > 0
        fprintf(1,'   Filtering Data\n');
        data = gaussianfilterdata(data,smoothSigma);
    end
        
    N = length(data);
    
    if nargin < 6 || isempty(maxNumber) || maxNumber <= 0
        maxNumber = N;
    end
    
    if s*threshold < minNoiseLevel
        s = 1;
        threshold = minNoiseLevel;
    end
        
    peakIdx = findWaveformPeaks(data,s,threshold,diffThreshold,maxNumber);
    peakIdx = peakIdx(peakIdx > diffThreshold & peakIdx < N - diffThreshold);
    peakIdx = peakIdx(data(peakIdx - diffThreshold) ~= 0 & data(peakIdx + diffThreshold) ~= 0);
    L = length(peakIdx);
    
    
    normalizedPeaks = zeros(L,2*diffThreshold+1);
    amplitudes = zeros(L,1);
    for i=1:L
        normalizedPeaks(i,:) = data(peakIdx(i)-diffThreshold:peakIdx(i)+diffThreshold);
        amplitudes(i) = sqrt(sum(normalizedPeaks(i,:).^2));
        normalizedPeaks(i,:) = normalizedPeaks(i,:) ./ amplitudes(i);
        if data(peakIdx(i)) < 0
            normalizedPeaks(i,:) = - normalizedPeaks(i,:);
        end
    end
    
    
    