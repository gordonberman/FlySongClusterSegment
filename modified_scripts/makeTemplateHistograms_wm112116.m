function makeTemplateHistograms_wm112116(templates,bins,means,ylimits,...
    isNoise,percentBelowNoiseThreshold)
    %makes template histogram plot
    
    %Inputs:
    
    %templates -> L x 1 cell array of templates (or groupings)
    %bins -> number of bins in each column (default = 50)
    %means -> lines to be plotted on top of histograms (the template mean
    %           if unspecified)
    %ylimits -> [ymin ymax] for plots (default = [-.5 .5])
    %isNoise (optional) -> isNoise from template creation for plot labels (added
    %11/11/16)
    %percentBelowNoiseThreshold (optional) -> percentBelowNoiseThreshold from
    %template creation for plot labels (added 11/17/16)

    N = length(templates);
    d = length(templates{1}(1,:));
    r = floor(d/2);
    snlabel = {'signal','noise'}; %added 11/11/16

    if nargin < 2 || isempty(bins)
        bins = 50;
    end

    if nargin < 3 || isempty(means) || length(means) ~= N
        means = cell(N,1);
        for i = 1:N
            means{i} = mean(templates{i});
        end
    end
    
    if nargin < 4 || isempty(ylimits)
        ylimits = [-.5 .5];
    end
    
    L = ceil(sqrt(N));
    M = ceil(N/L);
    qq = -r:r;
    if length(qq) > length(templates{1}(1,:))
        qq = qq(1:end-1);
    end
    
    for i=1:N
        
        subplot(M,L,i)
        
        xx = linspace(ylimits(1),ylimits(2),bins);
        Z = zeros(d,bins);
        for j=1:d
            Z(j,:) = hist(templates{i}(:,j),xx);
        end
        Z = bsxfun(@rdivide,Z,sum(Z,2))/(xx(2)-xx(1));
        
        pcolor(qq,xx,Z');
        shading flat
        hold on
        caxis([0 .02*diff(ylimits)])
        
        
        plot(qq,means{i},'k-','linewidth',2)
        if nargin==6 && ~isempty(isNoise) && ~isempty(percentBelowNoiseThreshold)
            title({['Template #' num2str(i) ', N = ' ...
                num2str(length(templates{i}(:,1)))], snlabel{isNoise(i)+1}, ...
                 ['(' num2str(percentBelowNoiseThreshold(i),3) ' < noiseThreshold)']}); %modified 11/17/16
        else
            title(['Template #' num2str(i) ', N = ' num2str(length(templates{i}(:,1)))]);
        end
        
    end
    
    addpath(genpath('./utilities'));
    load('saved_colormaps.mat','cc');
    colormap(cc)