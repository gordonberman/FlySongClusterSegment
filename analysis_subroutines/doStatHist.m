function doStatHist(pulsetrains, whichstat, options, outbase, Species)
    %Take the output of getSongParameters and summarize
    %Make histogram of whichstat, along with a text file
    %Keep format similar to previous versions
    %To do: make plotting optional
    %To do: optionally save 'alltGStat' (all values of whichstat)
    
    %Get the full options (including numIPIbins and IPI_sigma)
    options.setAll = false;
    options = makeParameterStructure(options);
    %Don't convert IPI_sigma because distances in pulsetrains are in ms
    %IPI_sigma = options.IPI_sigma * options.fs / 1000;
    
    %Sort by grouping, if desired
    %Why does this work but not pulsetrains.(whichstat)?
    %To do: Make templateGroups optional
    hastempgroups = false;
    if isfield(pulsetrains, 'templateGroups')
        hastempgroups = true;
        alltG = pulsetrains.('templateGroups');
        tG = unique(alltG);
    else
        tG = '';
        TemplateGroup = '';
        outbasetg = strjoin({outbase, whichstat},'_');
    end
    
    for g=1:length(tG)
        if hastempgroups
            TemplateGroup = tG(g);
            outbasetg = strjoin({outbase, whichstat, num2str(tG(g))},'_');
        end
        
        %Can't get indexing to work; instead trying a loop
        %ptIndex = find(alltG == tG(g));
        %Concatenate across all trains
        alltGstat = [];
        for pt = 1:length(pulsetrains)
            if hastempgroups
                if alltG(pt) == TemplateGroup
                    %Filter for length, if desired
                    if length(pulsetrains(pt).pulses) >= options.minpulse
                        alltGstat = [alltGstat pulsetrains(pt).(whichstat)];
                    end
                end
            else
                %Include all pulse trains above 'minpulse' threshold
                if length(pulsetrains(pt).pulses) > options.minpulse
                    alltGstat = [alltGstat pulsetrains(pt).(whichstat)];
                end
            end
        end
        
        %Exclude if sample size below minno
        if length(alltGstat) < options.minno
            [~,basename,~] = fileparts(outbasetg);
            fprintf(1,['Summary of %s not created for group %.0f due to ' ...
                'low sample size (%.0f) in: %s\n'], whichstat, ...
                TemplateGroup, length(alltGstat), basename);
            continue
        end
        %fprintf(1,'Total number of %s for group %.0f is %.0f\n',whichstat,...
        %    TemplateGroup,length(alltGstat));
        
        %Write all data
        save(strjoin({outbasetg, ['min' ...
            num2str(options.minpulse) 'pulse.mat']},'_'), 'alltGstat');
        dlmwrite(strjoin({outbasetg, ['min' ...
            num2str(options.minpulse) 'pulse.txt']},'_'),alltGstat,'delimiter','\n');
        
        %Make histogram plot (before or after smoothing?)
        %numIPIBins = floor(length(alltGstat)/10);
        %hist(alltGstat,numIPIBins);
        hist(alltGstat,100); %numIPIBins is too high for plotting
        savefig(strjoin({outbasetg, ['min' ...
            num2str(options.minpulse) 'pulse'], 'hist.fig'},'_'));
        close;
        
        %Run Gordon's mini-script to summarize (from createTemplates)
        %find IPI distribution through kernel density estimation
        [Y,X] = hist(alltGstat,options.numIPIBins);
        Y = normalizeHist(X,Y);
        IPI_sigma = options.IPI_sigma / (X(2) - X(1));
        Y = gaussianfilterdata(Y,IPI_sigma);   
        s = fit(X',Y','spline');

        %estimate mode of the IPI distribution
        StatMaxDensity = fminsearch(@(x) -s(x),X(argmax(Y)));
        MaxDensityPeakValue = s(StatMaxDensity);
        
        %get other summaries
        StatMean = mean(alltGstat);
        StatMedian = median(alltGstat);
        StatLQ = quantile(alltGstat,0.25);
        StatUQ = quantile(alltGstat,0.75);
        StatTot = sum(alltGstat);
        NoObs = length(alltGstat);
        
        %Write summary data
        if hastempgroups
            stattab = table({Species},TemplateGroup,StatMean,StatMedian,StatMaxDensity,...
                MaxDensityPeakValue,StatLQ,StatUQ,StatTot,NoObs);
        else
            stattab = table({Species},StatMean,StatMedian,StatMaxDensity,...
            MaxDensityPeakValue,StatLQ,StatUQ,StatTot,NoObs);
        end
        stattab.Properties.VariableNames{'Var1'} = 'Species';
        writetable(stattab,strjoin({outbasetg, ['min' ...
            num2str(options.minpulse) 'pulse'], 'summarystats.csv'},'_'));
    end
end