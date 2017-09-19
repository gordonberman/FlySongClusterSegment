function CreateAndAssignTemplatesOpts(trainfile,valfile,fscspath,minnoise,...
    ampthresh,noisepost,smoothnoise,numtemp,clusttemp,numsongtypes,...
    options,manualassign) 
    %Input:
    %trainfile (string): the training file, a 15+ minute .wav file with multiple
    %pulse trains per species
    %valfile (string): the validation file, a .wav file with peaks counted in each
    %subset (as made by randomizeAndCombineSignalNoiseWav.m)
    %fscspath (string): path to current version of FlySongClusterSegment
    %minnoise (numeric): min_noise_threshold, a minimum cutoff of the
    %smoothed data below which to mask as noise, can be set by converting
    %the minimum threshold in '_thresholds.csv' from createTrainingSet (default = -1
    %to ignore)
    %ampthresh (numeric): amplitude_threshold, percentage of peak amplitudes below the noise threshold 
    %to be called a "noise" template (default = .75)
    %smoothnoise (numeric): smoothingLength_noise, smoothing length for noise filter in milliseconds, 
    %should be about 1/2 the width of the smallest pulse (default = 4)
    %noisepost (numeric): noise_posterior_threshold (between 0 and 1, 1 most stringent, default = .5)
    %numtemp (numeric): number of templates to generate
    %clusttemp (logical/numeric): whether to cluster templates by song type
    %numsongtypes (numeric): number of expected song types
    %options (string): a .mat file containing a structure called 'options'
    %with default options (over-ridden by any other variables)
    %manualassign (string/numeric/logical): use manual assignment of
    %signal/noise?
    %Output:
    %fullfile(outdir1,[outbase1 '_TemplateHistograms.fig'])): template histograms
    %fullfile(outdir1,[outbase1 '_outputCreateTemplates.mat']): all output variables from createTemplates,
    %including outputData (containing all information about templates)
    %fullfile(outdir2,[outbase2 '_outputAssignTemplates.mat']): all output variables from assignDataToTemplates,
    %including peakIdxGroup (all called peaks by template number)
    %fullfile(outdir2,[outbase2 '_TemplatesAssigned.fig'])): a plot of all peak calls (assigned templates)
    %%Previously set: halfwidths (numeric): number of full widths at half maximum of IPI distribution
    %%to set for diffThreshold, should be related to smoothingLength_noise (default = 4)
    %%filtassign (1/0 or true/false): whether to high-pass filter the validation file before
    %%assigning templates

    %Set defaults
    if exist('options','var') %load options from options structure
        load(options);
        if isfield(options,'amplitude_threshold')
            ampthresh = options.amplitude_threshold;
        end
        if isfield(options,'smoothingLength_noise')
            smoothnoise = options.smoothingLength_noise;
        end
        if isfield(options,'noise_posterior_threshold')
            noisepost = options.noise_posterior_threshold;
        end
        if isfield(options,'min_noise_threshold')
            minnoise = options.min_noise_threshold;
        end
        if isfield(options,'k')
            numtemp = options.k;
        end
        if isfield(options,'humanLabel')
            manualassign = options.humanLabel;
        end
    end
    if nargin < 3 || isempty(fscspath)
        %either set a default path or assume that scripts are in the same directory
        fscspath = pwd;
    end
    if nargin < 4 || isempty(minnoise)
        minnoise = -1;
    end
    if nargin < 5 || isempty(ampthresh)
        ampthresh = 0.3;
    end
    if nargin < 6 || isempty(noisepost)
        noisepost = 0.5;
    end
    if nargin < 7 || isempty(smoothnoise)
        smoothnoise = 4;
    end
    if nargin < 8 || isempty(numtemp)
        numtemp = 12;
    end
    if nargin < 9 || isempty(clusttemp)
        clusttemp = false;
    end
    if nargin < 12 || isempty(manualassign)
        manualassign = false;
    end

    %Add path to FlySongClusterSegment
    addpath(genpath(fscspath));
    %Test to ensure scripts can be found
    try
        createTemplates([],[]);
    catch MEtest
        if strcmp(MEtest.identifier,'MATLAB:UndefinedFunction')
            error(['MATLAB cannot find your scripts.\nPlease specify' ...
                  ' the full path to FlySongClusterSegment in fscspath\n' ...
                  ' (the third argument to CreateAndAssignTemplatesOpts),\n'
                  ' or use addpath(genpath(yourpath)) to add this path.']);
        end
    end 

    %Convert any string variables to numeric
    if ischar(minnoise)
        minnoise = str2double(minnoise);
    end
    if ischar(ampthresh)
        ampthresh = str2double(ampthresh);
    end
    if ischar(smoothnoise)
        smoothnoise = str2double(smoothnoise);
    end
    if ischar(noisepost)
        noisepost = str2double(noisepost);
    end
    if ischar(numtemp)
        numtemp = str2double(numtemp);
    end
    if ischar(clusttemp)
        clusttemp = str2num(clusttemp);
    end
    if ischar(manualassign)
        manualassign = str2num(manualassign);
    end
    if exist('numsongtypes','var') == 1 && ischar(numsongtypes)
        numsongtypes = str2num(numsongtypes);
    end
    
    %Generate templates using the training file
    [outdir1,outbase1,~] = fileparts(trainfile);
    newoutdir1 = [outdir1 '/templates_' outbase1];
    if exist(newoutdir1,'dir') ~= 7
        try
            mkdir(newoutdir1);
        catch MEperm
            sprintf('Unable to create subdirectory; storing output in %s.\n',...
                outdir1);
            newoutdir1 = outdir1;
        end
    end
    outbase1 = [outbase1 '_at' num2str(ampthresh) '_sn' num2str(smoothnoise) ...
        '_np' num2str(noisepost) '_mn' num2str(minnoise) '_nt' num2str(numtemp)];
    try
        [d,fs] = audioread(trainfile);
    catch MEnotrain
        fprintf(1,'%s missing\n',trainfile);
        rethrow(MEnotrain);
    end

    options.fs = fs;
    options.amplitude_threshold = ampthresh; %default .75
    options.smoothingLength_noise = smoothnoise; %default 4
    options.noise_posterior_threshold = noisepost; %default 0.5; vary this for ROC curve
    options.min_noise_threshold = minnoise; %default -1; use sqrt(10^(minimum
    %of thresholds from 'createTrainingSet')
    options.k = numtemp; %number of templates to generate; increase if 
    %signal peaks get called with noise templates
    if isfield(options,'min_seperation')
        if options.min_seperation > 2
            options.min_seperation = 1; %reduce minimum distance between pulses
            %(should be unnecessary when Gordon updates FSCS)
        end
    end
    options.humanLabel = manualassign; %default true (?)
    %new: try Gordon's t-SNE analysis... this takes a while!
    %options.run_tsne = true;
    [outputData,allPeakIdx,allNormalizedPeaks,peakAmplitudes,isNoise,allScores,...
        options] = createTemplates_wm112116(d,options);
    %[outputData,allPeakIdx,allNormalizedPeaks,peakAmplitudes,isNoise,allScores,...
    %    options] = createTemplates(d,options);
    %tunning t-SNE creates 2 plots; save both
    if options.run_tsne
        savefig(fullfile(newoutdir1,[outbase1 '_tSNEPlot.fig']));
        close;
    end
    savefig(fullfile(newoutdir1,[outbase1 '_TemplateHistograms.fig']));
    close;
    save(fullfile(newoutdir1,[outbase1 '_outputCreateTemplates.mat']), 'outputData',...
        'allPeakIdx','allNormalizedPeaks','peakAmplitudes','isNoise','allScores',...
        'options');
    
    %Now assign the validation data to templates
    %Make the output filenames match those from AssignTemplatesAndGetFPFN
    %(i.e., preserve info about the filename from training/templates).
    tempsuff = ['temp' outbase1];
    valsuff = ['valnp' num2str(noisepost)];
    [outdir2,outbase2,~] = fileparts(valfile);
    [dfilt2,fs2] = audioread(valfile);
    if length(dfilt2) > 600*fs2
        %reduce length of valfile to 10 minutes
        dfilt2 = dfilt2(1:(600*fs2));
    end
    newoutdir2 = [outdir2 '/output_' outbase2];
    if exist(newoutdir2,'dir') ~= 7
        try
            mkdir(newoutdir2);
        catch MEperm2
            sprintf('Unable to create subdirectory; storing output in %s.\n',...
                outdir2);
            newoutdir2 = outdir2;
        end 
    end
    outbase2 = strjoin({outbase2, tempsuff, valsuff},'_');
    options.fs = fs2;
    [groupings,peakIdxGroup,likes,allPeakIdx,allNormalizedPeaks,noiseThreshold] = ...
                     assignDataToTemplates(dfilt2,outputData,options);
    save(fullfile(newoutdir2,[outbase2 '_outputAssignTemplates.mat']), 'groupings',...
        'peakIdxGroup','likes','allPeakIdx','allNormalizedPeaks','noiseThreshold');
    %Plot the templates IF there aren't a crazy number of peaks
    if (length(allPeakIdx) < 5700)
        makePeakPlot(dfilt2,peakIdxGroup,1:length(peakIdxGroup));
        savefig(fullfile(newoutdir2,[outbase2 '_TemplatesAssigned.fig']));
    else
        fprintf(1, '    Too many peaks to plot!\n');
    end
    %Plot only signal templates (again, if not too many)
    signalPeakIdx = getSignalPeakIdx(peakIdxGroup, isNoise);
    if length(signalPeakIdx) < 5700
        makePeakPlot(dfilt2,peakIdxGroup,find(1-isNoise));
        savefig(fullfile(newoutdir2,[outbase2 '_SignalTemplatesAssigned.fig']));
    else
        fprintf(1, '    Too many peaks to plot!\n');
    end
    
    %Do both clustered and un-clustered analysis
    if clusttemp
        outbase1 = [outbase1 '_clustered'];
        outbase2 = [outbase2 '_clustered'];
        if sum(isNoise) > 0
            numsongtypes = numsongtypes + 1;
        end
        outputData = clusterTemplates_wm112116(outputData,numsongtypes,...
            options);
        %also save the clustered outputData!
        save(fullfile(newoutdir1,[outbase1 '_outputCreateTemplates.mat']), 'outputData',...
        'allPeakIdx','allNormalizedPeaks','peakAmplitudes','isNoise','allScores',...
        'options');
        for c=1:numsongtypes
            savefig(fullfile(newoutdir1,[outbase1 '_TemplateHistograms_group' ...
                num2str(numsongtypes-c+1) '.fig']));
            close;
        end
        [groupings,peakIdxGroup,likes,allPeakIdx,allNormalizedPeaks,noiseThreshold] = ...
                         assignDataToTemplates(dfilt2,outputData,options);
        save(fullfile(newoutdir2,[outbase2 '_outputAssignTemplates.mat']),...
            'groupings','peakIdxGroup','likes','allPeakIdx',...
            'allNormalizedPeaks','noiseThreshold');
        %Plot the templates IF there aren't a crazy number of peaks
        if (length(allPeakIdx) < 5700)
            makePeakPlot(dfilt2,peakIdxGroup,1:length(peakIdxGroup));
            savefig(fullfile(newoutdir2,[outbase2 '_TemplatesAssigned.fig']));
        else
            fprintf(1, '    Too many peaks to plot!\n');
        end
        %Plot only signal templates (again, if not too many)
        signalPeakIdx = getSignalPeakIdx(peakIdxGroup, outputData.isNoiseTemplateGrouping);
        if length(signalPeakIdx) < 5700
            makePeakPlot(dfilt2,peakIdxGroup,find(1-outputData.isNoiseTemplateGrouping));
            savefig(fullfile(newoutdir2,[outbase2 '_SignalTemplatesAssigned.fig']));
        else
            fprintf(1, '    Too many peaks to plot!\n');
        end
    end
end
