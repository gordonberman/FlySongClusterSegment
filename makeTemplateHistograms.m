function makeTemplateHistograms(templates,bins,means)
    %makes template histogram plot
    
    %Inputs:
    
    %templates -> L x 1 cell array of templates (or groupings)
    %bins -> number of bins in each column (default = 50)
    %means -> lines to be plotted on top of histograms (the template mean
    %           if unspecified)

    N = length(templates);
    d = length(templates{1}(1,:));

    if nargin < 2 || isempty(bins)
        bins = 50;
    end

    if nargin < 3 || isempty(means) || length(means) ~= N
        means = cell(N,1);
        for i = 1:N
            means{i} = mean(templates{i});
        end
    end
    
    
    L = ceil(sqrt(N));
    M = ceil(N/L);
    for i=1:N
        
        subplot(M,L,i)
        
        xx = linspace(-.5,.5,bins);
        Z = zeros(d,bins);
        for j=1:d
            Z(j,:) = hist(templates{i}(:,j),xx);
        end
        pcolor(1:d,xx,Z');
        shading flat
        hold on
        
        plot(1:d,means{i},'k-','linewidth',2)
        title(['Template #' num2str(i) ', N = ' num2str(length(templates{i}(:,1)))]);
        
    end