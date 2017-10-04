function [out,mask,threshold,smoothedPeakLocations] = filterDataAmplitudes(...
    data,sigma,minLength,maxNumGaussians,replicates,maxNumPoints,minSeperation,...
    threshold,min_threshold,noise_posterior_threshold)

    if nargin < 2 || isempty(sigma)
        sigma = 40;
    end
    
    if nargin < 3 || isempty(minLength)
        minLength = 20;
    end

    if nargin < 4 || isempty(maxNumGaussians)
        maxNumGaussians = 2;
    end
    
    if nargin < 5 || isempty(replicates)
        replicates = 3;
    end
    
    if nargin < 6 || isempty(maxNumPoints)
        maxNumPoints = 10000;
    end
    
    if nargin < 7 || isempty(minSeperation)
        minSeperation = 100;
    end
    
    if nargin < 9 || isempty(min_threshold)
        min_threshold = -9e99;
    end
    
    if nargin < 10 || isempty(noise_posterior_threshold)
        noise_posterior_threshold = .5;
    end
    
    s = size(data);
    if s(2) > s(1)
        data = data';
    end
    
    %log of gaussian filter smoothed squared signal
    y = log10(gaussianfilterdata(data.^2,sigma));
    
    %find threshold and initial mask    
    if nargin < 8 || isempty(threshold) 
    
        %gaussian mixture model for log(smoothed signal)
        obj = findBestGMM_AIC(y,maxNumGaussians,replicates,maxNumPoints);
        
        %identify max peak and min peak
        idx = argmax(obj.mu);
        minIdx = argmin(obj.mu);
        posts = posterior(obj,y);
        posts = posts(:,idx);
        
        %find mask and threshold
        try
            mask = (posts >= noise_posterior_threshold | y > obj.mu(idx)) & y > obj.mu(minIdx);
        catch MEdimineq
            size(posts)
            size(noise_posterior_threshold)
            noise_posterior_threshold
            rethrow(MEdimineq)
        end
        threshold = min(y(mask));
        
        if threshold < min_threshold
            threshold = min_threshold;
            mask = y > min_threshold;
        end
        
    else
        
        %find mask
        mask = y > threshold;
        
    end
    
    
    
    
    %zero all connected components smaller than minLength
    CC = bwconncomp(mask);
    lengths = returnCellLengths(CC.PixelIdxList);
    idx = find(lengths < minLength);
    for i=1:length(idx)
        mask(CC.PixelIdxList{idx(i)}) = false;
    end
    
    
    %make masked data set
    out = data.*double(mask);
    
    
    %find maxima of masked data set
    smoothedPeakLocations = find(imregionalmax(y.*double(mask)) & mask);
    
    
    %eliminate all peaks spaced less than diffThreshold
    if minSeperation ~= -1
        
        d = diff(smoothedPeakLocations);
        a = d < minSeperation;
        while sum(a) > 0
            
            idx = argmin(d);
            idx = idx(1);
            vals = [y(smoothedPeakLocations(idx)) y(smoothedPeakLocations(idx+1))];
            idx2 = argmin(vals);
            smoothedPeakLocations = setdiff(smoothedPeakLocations,smoothedPeakLocations(idx + idx2 - 1));
            
            d = diff(smoothedPeakLocations);
            a = d < minSeperation;
            
        end
    end
    
    
%     for i=1:length(smoothedPeakLocations)
%         dIdx = argmax(abs(data(smoothedPeakLocations(i) + (-sigma:sigma)))) - sigma;
%         smoothedPeakLocations(i) = smoothedPeakLocations(i) + dIdx;
%     end
    
    
%     if minSeperation ~= -1
%         
%         d = diff(smoothedPeakLocations);
%         a = d < minSeperation;
%         while sum(a) > 0
%             
%             if sum(a) > 0
%                 CC = bwconncomp(a);
%                 for j=1:CC.NumObjects
%                     b = [CC.PixelIdxList{j};CC.PixelIdxList{j}(end)+1];
%                     c = smoothedPeakLocations(b);
%                     minLocation = argmin(y(c));
%                     smoothedPeakLocations(b(minLocation)) = -1;
%                 end
%             end
%             
%             smoothedPeakLocations = smoothedPeakLocations(smoothedPeakLocations > 0);
%             d = diff(smoothedPeakLocations);
%             a = d < minSeperation;
%             
%         end
%         
%     end
    
    
    
