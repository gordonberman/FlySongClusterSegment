function [normalizedPeaks,peakIdx,data,amplitudes] = ...
    findNormalizedPeaks(data,s,threshold,diffThreshold,highpass,maxNumber,minNoiseLevel,samplingFreq)

    
    if nargin < 3 || isempty(threshold)
        threshold = 5;
    end
    
    if nargin < 4 || isempty(diffThreshold)
        diffThreshold = 150;
    end
    
    if nargin < 5 || isempty(highpass)
        highpass = 200;
    end
    
    if nargin < 7 || isempty(minNoiseLevel)
        minNoiseLevel = 0;
    end
    
    if nargin < 8 || isempty(samplingFreq)
        samplingFreq = 1e4;
    end
    
    %     if smoothSigma > 0
    %         fprintf(1,'   Filtering Data\n');
    %         data = gaussianfilterdata(data,smoothSigma);
    %     end
    
    %     if highpass > 0
    %         fprintf(1,'   Filtering Data\n');
    %         poles = 4;
    %         [b,a]=butter(poles,highpass./(samplingFreq/2),'low');
    %         data(:) = filtfilt(b,a,data);
    %     end
    
    [data,~,~] = filterDataAmplitudes(data,10,4,4,5,10000);
    
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
    
    
    