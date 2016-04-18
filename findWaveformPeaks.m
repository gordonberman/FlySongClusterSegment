function peakIdx = findWaveformPeaks(data,s,threshold,diffThreshold,keepNegatives,maxNumPeaks)

    if nargin < 3 || isempty(threshold)
        threshold = 3;
    end
    
    if nargin < 4 || isempty(diffThreshold)
        diffThreshold = 200;
    end
    
    if nargin < 5 || isempty(keepNegatives)
        keepNegatives = true;
    end
    
    N = length(data);
    
    if nargin < 6 || isempty(maxNumPeaks) || maxNumPeaks <= 0
        maxNumPeaks = N;
    end

    if keepNegatives
        peakIdx = find((regionalmax(data')|regionalmin(data')) & abs(data') > threshold*s);
    else
        peakIdx = find(regionalmax(data') & abs(data') > threshold*s);
    end
    vals = abs(data(peakIdx));
    
    reg_maxes = false(size(peakIdx));
    count = 0;
    for i=1:length(peakIdx)
        if vals(i) == max(vals(abs(peakIdx - peakIdx(i)) <= diffThreshold))
            reg_maxes(i) = true;
            count = count + 1;
            
            if count >= maxNumPeaks
                break;
            end
        end
    end
        
    peakIdx = peakIdx(reg_maxes);