function [obj,AICs] = findBestGMM_AIC(data,maxPeaks,replicates,maxNum)


    if nargin < 3 || isempty(replicates)
        replicates = 1;
    end

    
    N = length(data(:,1));
    if nargin < 4 || isempty(maxNum)
        maxNum = N;
    end
    
    if maxNum < N
        idx = randperm(N,maxNum);
        data = data(idx,:);
    end
    
    
    AICs = zeros(maxPeaks,1);
    objs = cell(maxPeaks,1);
    for i=1:maxPeaks
        %objs{i} = gmdistribution.fit(data,i,'Options',options,'Replicates',replicates,'Regularize',1e-30);
        objs{i} = gmixPlot(data,i,[],[],true,[],[],[],replicates);
        x = objs{i}.AIC;
        k = objs{i}.NDimensions;
        numParameters = i-1 + k*i + i*k*(k+1)/2;
        AICs(i) = x + 2*(numParameters*(numParameters+1))/(length(data(:,1)) - numParameters - 1);
    end
    
    minIdx = argmin(AICs);
    
    obj = objs{minIdx};