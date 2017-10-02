function [xMax,p] = findExtremumParabola(x,y,numPoints)

    if nargin < 3 || isempty(numPoints)
        numPoints = 1;
    end
    
    N = length(y);
    [~,maxIdx] = max(y);
    
    idx = (-numPoints:numPoints) + maxIdx;
    if min(idx) < 1
        idx = idx - min(idx) + 1;
    end
    if max(idx) > N
        idx = idx - (max(idx) - N);
    end
    
    p = polyfit(x(idx),y(idx),2);
    
    xMax = -.5*p(2)/p(1);
    