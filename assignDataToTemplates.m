function [groupings,peakIdxGroup,likes,allPeakIdx,allNormalizedPeaks,coeffs,projStds,baselines] = ...
                     assignDataToTemplates(data,templates,isNoise,options)

                 
    %Inputs:
                 
    %data -> 1-D array of data points
    %templates -> L x 1 cell array containing M_i x d arrays of template peaks
    %options -> self-explanatory

    
    %Outputs:
    
    %groupings -> L x 1 cell array, each cell containing an N_i x d array 
    %               of grouped peaks
    %peakIdxGroup -> L x 1 cell array of index values for each 
    %               corresponding group
    %likes -> N x L array of the log-likelihood values that each peak 
    %           belongs to the corresponding model
    %allPeakIdx -> all found peak locations
    %allNormalizedPeaks -> N x d array containing all found peaks 
    %                       corresponding to the locations in allPeakIdx
    %coeffs -> L x 1 cell array of model PCA bases
    %projStds -> L x 1 cell array of model projection standard deviations
    
    
    addpath(genpath('./chronux'))
    
    L = length(templates);
    d = length(templates{1}(1,:));
       
    if nargin < 4 || isempty(options)
        options.setAll = true;
    else
        options.setAll = false;
    end
    options = makeDefaultOptions(options);
    options.diffThreshold = (d-1) / 2;
    
    if options.noiseLevel <= 0
        fprintf(1,'   Finding Noise Level\n');
        if numel(data) > 1e6
            shortdata = data(1:1e6);
        else
            shortdata = data;
        end
        if options.smoothSigma > 0
            shortdata = gaussianfilterdata(shortdata,options.smoothSigma);
        end
        [ssf] = sinesongfinder(shortdata,options.fs,12,20,.1,.01,.05,[0 1000]);
        %returns ssf, which is structure containing the following fields
        noise = findnoise(ssf,3,80,1000);
        options.noiseLevel = noise.sigma;
    end
    
    %find all potential peaks in the new data set
    fprintf(1,'   Finding Preliminary Peak Locations\n');
    [normalizedPeaks,peakIdx,~] = ...
        findNormalizedPeaks(data,options.noiseLevel,options.sigmaThreshold,...
        options.diffThreshold,options.smoothSigma);
    
    N = length(normalizedPeaks(:,1));
    allPeakIdx = peakIdx;
    allNormalizedPeaks = normalizedPeaks;
    
    
    
    fprintf(1,'   Finding Template Bases and Projections\n');
    coeffs = cell(L,1);
    projStds = cell(L,1);
    means = cell(L,1);
    L_templates = zeros(L,1);
    
    for i=1:L
        
        fprintf(1,'      Template #%2i\n',i);
        L_templates(i) = length(templates{i}(:,1));
        
        %find Data Set Mean
        means{i} = mean(templates{i});
    
        %perform PCA on set of normalized peaks
        [coeffs{i},scores,~] = princomp(templates{i});
        
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
    
    %find baseline noise levels
    fprintf(1,'   Calculating Baseline Noise Levels\n');
    if options.use_likelihood_threshold
        baselines = repmat(options.baseline_threshold,sum(~isNoise),1);
    else
        [baselines,~] = findTemplateBaselines(templates,coeffs,means,projStds,isNoise,options.baseline_quantile);
    end
        
    %find likelihood values
    fprintf(1,'   Finding Likelihoods\n');
    likes = zeros(N,L);
    
    for i=1:L
        
        fprintf(1,'      Template #%2i\n',i);
        
        projections = (normalizedPeaks - repmat(means{i},N,1))*coeffs{i};
                
        likes(:,i) = findDataSetLikelihoods(projections,zeros(size(means{i})),projStds{i});
    end
    
    clear projections 
    
    
    [~,idx] = max(likes,[],2);
    
    %find assignments that are sub-baseline and assign them to to either
    %the next super-baseline signal class or the largest noise class
    fprintf(1,'   Checking Peaks Against Baseline\n');
    notNoise = find(~isNoise);
    %noiseVals = find(isNoise);
    signalLikes = likes(:,~isNoise);
    belowThreshold = signalLikes < repmat(baselines',length(idx),1);
    %[sortVals,sortIdx] = sort(likes,2,'descend');
    
    for q=1:length(notNoise)
        
        i = notNoise(q);
        
        idx2 = find(idx==i & belowThreshold(:,i));
        if ~isempty(idx2)
            for j=1:length(idx2)
                
                qq = likes(idx2(j),notNoise);
                qq(likes(idx2(j),notNoise) - baselines' < 0) = -10;
                [maxVal,maxIdx] = max(qq);
                if maxVal < 0
                    idx(idx2(j)) = L + 1;
                else
                    idx(idx2(j)) = notNoise(maxIdx);
                end
                
             
            end
        end
    end
    

    if max(idx) == L+1
        L = L+1;
    end
    
    groupings = cell(L,1);
    peakIdxGroup = cell(L,1);
    for i=1:L
        q = idx == i;
        groupings{i} = allNormalizedPeaks(q,:);
        peakIdxGroup{i} = allPeakIdx(q);
    end
    
    
    
    
    