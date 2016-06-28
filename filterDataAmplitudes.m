function [out,mask,obj,smoothedPeakLocations] = filterDataAmplitudes(data,sigma,minLength,maxNumGaussians,replicates,maxNumPoints)

    if nargin < 2 || isempty(sigma)
        sigma = 10;
    end

    y = gaussianfilterdata(abs(data),sigma);
    
    obj = findBestGMM_AIC(log10(y),maxNumGaussians,replicates,maxNumPoints);
    
    idx = argmax(obj.mu);
    posts = posterior(obj,log10(y));
    posts = posts(:,idx);
    
    mask = posts >= .5;
    CC = bwconncomp(mask);
    lengths = returnCellLengths(CC.PixelIdxList);
    idx = find(lengths < minLength);
    for i=1:length(idx)
        mask(CC.PixelIdxList{idx(i)}) = false;
    end
    
    out = data.*double(mask);
    
    smoothedPeakLocations = find(imregionalmax(y.*double(mask)));