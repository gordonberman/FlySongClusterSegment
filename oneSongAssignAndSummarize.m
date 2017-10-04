function oneSongAssignAndSummarize(songfile, tempfile, Species, ...
    outdir, summmaxIPI, minno, minpulse, plotassign, noisepost) 
    %To do: save plot of assigned templates
    %To do: split outputs into subdirectories
    %To do: make some outputs and plots optional
    addpath(genpath(pwd));

    %Set defaults for options
    if nargin < 5 || isempty(summmaxIPI)
        summmaxIPI = 200; 
    end
    if nargin < 6 || isempty(minno)
        minno = 100; 
    end
    if nargin < 7 || isempty(minpulse)
        minpulse = 5; 
    end
    if nargin < 8 || isempty(plotassign)
        plotassign = false;
    end
    if nargin < 9 || isempty(noisepost)
        noisepost = 0.5; %Which should take priority, the options or the input?
                         %This should be different from the options from
                         %generating templates
    end
    
    %Create outdir, if it doesn't exist
    if exist(outdir, 'dir') ~= 7
        try
            mkdir(outdir);
        catch MEperm
            fprintf(1,'Cannot make directory %s; exiting\n',outdir);
            rethrow(MEperm);
        end
    end
    
    %Make basename for output from song, template, and options names
    [~,songbase,~] = fileparts(songfile);
    [~,tempbase,~] = fileparts(tempfile);
    outbase = strjoin({songbase, tempbase},'_');
    optsdir = fullfile(outdir, 'optionsused');
    if exist(optsdir,'dir') ~= 7
        try
            mkdir(optsdir);
        catch MEperm2
            fprintf(1,['Cannot make options directory\n\(%s\);\nSaving options' ...
                'files to output directory:\n%s\n'],optsdir,outdir);
            optsdir = outdir;
        end
    end
    %optsfile = fullfile(outdir, [outbase '_options.mat']); %save for records
    optsfile = fullfile(optsdir, [outbase '_options.mat']); %save for records

    %Load templates, including options
    sprintf('Loading templates: %s\n',tempfile);
    load(tempfile, 'outputData', 'isNoise', 'options');
    
    %Set new options AFTER loading options from templates
    options.summmaxIPI = summmaxIPI; %maximum distance between pulses in ms 
    options.minno = minno; %minimum total number of observations (e.g. pulses) to include a recording
    options.minpulse = minpulse; %minimum number of pulses in a pulse train
    options.numIPIBins = 10000; %for filterDataAmplitudes
    options.IPI_sigma = 1; %selected by eye from plot for IPI
    options.noise_posterior_threshold = noisepost; %can be different from threshold 
    %set when creating templates
    %separately optimize for PTL, PPB, and CF
    save(optsfile,'options');
    
    %Read in songfile
    sprintf('Loading song file %s\n',songfile);
    [song, fs] = audioread(songfile);
    if fs ~= options.fs
        song = fixSamplingFrequency(song,fs,options.fs);
    end
    
    %Assign data to templates
    if exist(fullfile(outdir,[outbase '_outputAssignTemplates.mat']),...
            'file') ~= 2
        sprintf('Assigning data to templates\n');
        [groupings,peakIdxGroup,likes,allPeakIdx,allNormalizedPeaks,...
            noiseThreshold,freqIdxGroup] = assignDataToTemplates(song,outputData,options);
        save(fullfile(outdir,[outbase '_outputAssignTemplates.mat']), 'groupings',...
            'peakIdxGroup','likes','allPeakIdx','allNormalizedPeaks','noiseThreshold',...
            'freqIdxGroup');
        if plotassign
            %make and save plot of assignment
            if exist(fullfile(outdir,[Species 'Assignments'])) ~= 7
                mkdir(fullfile(outdir,[Species 'Assignments']));
            end
            if isfield(outputData,'templateGroupings')
                makePeakPlot(song,peakIdxGroup,find(1-outputData.isNoiseTemplateGrouping));
            else
                makePeakPlot(song,peakIdxGroup,find(1-outputData.isNoise));
            end
            savefig(fullfile(outdir,[Species 'Assignments'], ...
                [outbase '_SignalTemplatesAssigned.fig']));
            close;
        end
    else
        load(fullfile(outdir,[outbase '_outputAssignTemplates.mat']),...
            'peakIdxGroup','freqIdxGroup');
    end
                 
    %Get pulse trains
    if exist(fullfile(outdir,[outbase '_pulseTrains.mat']),'file') ~= 2
        if isfield(outputData,'templateGroupings')
            %pulsetrains = getSongParameters(song, peakIdxGroup, isNoise, ...
            %    outputData.templateGroupings, options)
            [pulsetrains, options] = getSongParameters(songbase, ...
                peakIdxGroup, outputData.isNoiseTemplateGrouping, true, ...
                options, freqIdxGroup);
        else
            [pulsetrains, options] = getSongParameters(songbase, ...
                peakIdxGroup, isNoise, false, options, freqIdxGroup);
        end
        save(fullfile(outdir,[outbase '_pulseTrains.mat']),'pulsetrains');
    else
        load(fullfile(outdir,[outbase '_pulseTrains.mat']),'pulsetrains');
    end
    
    %Summarize and plot histogram (make plotting optional in doStatHist)
    %Also allow unclustered data
    if isfield(outputData,'templateGroupings')
        outsuff = strjoin({num2str(min(outputData.templateGroupings)), ['min' ...
            num2str(minpulse)], 'hist.fig'},'_');
    else
        outsuff = strjoin({['min' num2str(minpulse)], 'hist.fig'},'_');
    end
    if exist(strjoin({outbase, 'ipi', outsuff},'_'),'file') ~= 2
        doStatHist(pulsetrains, 'ipi', options, fullfile(outdir,outbase), Species);
    end
    if exist(strjoin({outbase, 'numPulses', outsuff},'_'),'file') ~= 2
        doStatHist(pulsetrains, 'numPulses', options, fullfile(outdir,outbase), Species);
    end
    if exist(strjoin({outbase, 'trainLengths', outsuff},'_'),'file') ~= 2
        doStatHist(pulsetrains, 'trainLengths', options, fullfile(outdir,outbase), Species);
    end
    if exist(strjoin({outbase, 'cf', outsuff},'_'),'file') ~= 2
        doStatHist(pulsetrains, 'cf', options, fullfile(outdir,outbase), Species);
    end
end
