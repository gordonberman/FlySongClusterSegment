function [newTemplates,isNoise,splitted] = selectTemplates(templates,plotsOn)
    %interface for selecting templates from data

    if nargin < 2 || isempty(plotsOn)
        plotsOn = true;
    end

    splitted = false;
    
    maxTemplates = 100;
    histogramBins = 50;
    samplingRate = 1e-4;
    
    N = length(templates);
    d = length(templates{1}(1,:));
    
    figure
    makeTemplateHistograms(templates,histogramBins);
   
    
    newTemplates = cell(maxTemplates,1);
    noiseTemplates = cell(maxTemplates,1);
    isNoise = true(maxTemplates,1);
    count = 1;
    countNoise = 1;
    for i=1:N
        
        figure(1101)
        
        xx = linspace(-.5,.5,histogramBins);
        Z = zeros(d,histogramBins);
        for j=1:d
            Z(j,:) = hist(templates{i}(:,j),xx);
        end
        pcolor((1:d).*samplingRate,xx,Z');
        shading flat
        hold on
        plot((1:d).*samplingRate,mean(templates{i}),'k-','linewidth',2)
        title(['Template #' num2str(i) ', N = ' num2str(length(templates{i}(:,1))) ...
            ', Signal (s), Noise (n), or Split (p)? ']);
        
        template_action = '1';
        while max(template_action ~= 's') && max(template_action ~= 'n') && max(template_action ~= 'p')
            [~,~,template_action] = ginput(1);
        end
        
             
        if template_action == 'n'
            noiseTemplates{countNoise} = templates{i};
            countNoise = countNoise + 1;
        end
        
        
        if template_action == 's'
            newTemplates{count} = templates{i};
            count = count + 1;
        end
        
        
        if template_action == 'p'
            idx = kmeans(templates{i},2);
            newTemplates{count} = templates{i}(idx==1,:);
            count = count + 1;
            newTemplates{count} = templates{i}(idx==2,:);
            count = count + 1;
            splitted = true;
        end
        
        figure(1101)
        hold off
        
    end
    
    close 1101
    
    newTemplates = newTemplates(1:count-1);
    isNoise(1:count-1) = false;
    
    noiseTemplates = noiseTemplates(1:countNoise-1);
    newTemplates = [newTemplates; noiseTemplates];
    isNoise = isNoise(1:length(newTemplates));
    
    
    if plotsOn
        figure
        
        makeTemplateHistograms(newTemplates,histogramBins);
        
    end
    