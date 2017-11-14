function doStatHist(pulsetrains, whichstat, options, outbase, Species,...
    doplots)
    %Testing 01/14/17: try getting the max density nearest the median
    %Take the output of getSongParameters and summarize
    %Make histogram of whichstat, along with a text file
    %Keep format similar to previous versions
    %To do: make plotting optional
    %To do: optionally save 'alltGStat' (all values of whichstat)
    if nargin < 6 || isempty(doplots)
        doplots = true;
    end
    
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
        %For cf, exclude pulses exactly at repetitive frequencies
        if strcmp(whichstat,'cf')
            if isfield(options,'repfreq')
                %repfreq should represent the base frequency
                %create an array of multiples of repfreq
                freqtoexc = [options.repfreq round(options.repfreq*2,2) ...
                    round(options.repfreq*3,2) round(options.repfreq*4,2) ...
                    round(options.repfreq*5,2)];
                olength = length(alltGstat);
                %exclude these values from alltGstat
                alltGstat = setdiff(alltGstat, freqtoexc);
                if length(alltGstat) < olength
                    fprintf(1, ['Removed %d carrier frequency values ' ...
                        'as multiples of repfreq (%d)\n'], ...
                        olength-length(alltGstat), options.repfreq);
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
        if doplots
            hist(alltGstat,100); %numIPIBins is too high for plotting
            savefig(strjoin({outbasetg, ['min' ...
                num2str(options.minpulse) 'pulse'], 'hist.fig'},'_'));
            close;
        end
        %Run Gordon's mini-script to summarize (from createTemplates)
        %find IPI distribution through kernel density estimation
        [Y,X] = hist(alltGstat,options.numIPIBins);
        if strcmp(whichstat,'cf')
            Y = medfilt1(Y,3);  %filter to eliminate regular noise
        end
        Y = normalizeHist(X,Y);
        IPI_sigma = options.IPI_sigma / (X(2) - X(1));
        Y = gaussianfilterdata(Y,IPI_sigma); 
        finidx = isfinite(X) & isfinite(Y); %to remove issue with NaNs
        %also plot output?
        if doplots && sum(finidx) > 2
            plot(X,Y);
            savefig(strjoin({outbasetg, ['min' ...
                num2str(options.minpulse) 'pulse'], 'smoothed.fig'},'_'));
            close;
        end
        if sum(finidx) > 2
            try
                s = fit(X(finidx)',Y(finidx)','spline');
            catch MEfit
                length(finidx)
                sum(~isnan(X))
                sum(~isnan(Y))
                rethrow(MEfit);
            end
            %estimate mode of the IPI distribution
            StatMaxDensity = fminsearch(@(x) -s(x),X(argmax(Y)));
            MaxDensityPeakValue = s(StatMaxDensity);
            %also estimate local maximum nearest median
            %(may sometimes be more accurate than starting at mode of smoothed distribution)
            StatMaxDensityNearMedian = fminsearch(@(x) -s(x),...
                median(alltGstat,'omitnan'));
            MaxDensityNearMedianPeakValue = s(StatMaxDensityNearMedian);
            StatMaxDensityNearSmoothedMedian = fminsearch(@(x) -s(x),...
                median(X,'omitnan'));
            MaxDensityNearSmoothedMedianPeakValue = s(StatMaxDensityNearSmoothedMedian);
        else
            if hastempgroups
            fprintf(1,['Insufficient smoothed data for estimating ' ...
                'mode for group %i for:\n%s\n'], TemplateGroup, outbase);
            else
                fprintf(1,['Insufficient smoothed data for estimating ' ...
                'mode for:\n%s\n'], outbase);
            end
            StatMaxDensity = NaN;
            MaxDensityPeakValue = NaN;
            StatMaxDensityNearMedian = NaN;
            MaxDensityNearMedianPeakValue = NaN;
            StatMaxDensityNearSmoothedMedian = NaN;
            MaxDensityNearSmoothedMedianPeakValue = NaN;  
        end
        
        %get other summaries
        StatMean = mean(alltGstat,'omitnan');
        StatMedian = median(alltGstat,'omitnan');
        StatLQ = quantile(alltGstat,0.25);
        StatUQ = quantile(alltGstat,0.75);
        StatTot = sum(alltGstat,'omitnan');
        NoObs = length(alltGstat);
        
        %Write summary data
        if hastempgroups
            stattab = table({Species},TemplateGroup,StatMean,StatMedian,StatMaxDensity,...
                MaxDensityPeakValue,StatMaxDensityNearMedian,...
                MaxDensityNearMedianPeakValue,StatMaxDensityNearSmoothedMedian,...
                MaxDensityNearSmoothedMedianPeakValue,StatLQ,StatUQ,StatTot,NoObs);
        else
            stattab = table({Species},StatMean,StatMedian,StatMaxDensity,...
            MaxDensityPeakValue,StatMaxDensityNearMedian,...
            MaxDensityNearMedianPeakValue,StatMaxDensityNearSmoothedMedian,...
            MaxDensityNearSmoothedMedianPeakValue,StatLQ,StatUQ,StatTot,NoObs);
        end
        stattab.Properties.VariableNames{'Var1'} = 'Species';
        writetable(stattab,strjoin({outbasetg, ['min' ...
            num2str(options.minpulse) 'pulse'], 'summarystats.csv'},'_'));
    end
end
