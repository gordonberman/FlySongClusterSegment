function SeparateSongsAssignTemplatesAndGetFPFN(tempfile,valsignaldir,...
    valnoisedir,categorytimes,varval,filelist,varname,varpref)
    %New 12/9/16: allow tweaking any input variable
    %varval: value, varname: name, varpref: prefix for filename
    %New 12/14/16: Assign templates separately to each file (because noise
    %levels vary across files)
    %Convert any string variables to numeric
    if ischar(varval)
        varval = str2double(varval);
    end
    if nargin < 8
        %make compatible with previous version
        varname = 'noise_posterior_threshold';
        varpref = 'valnp';
    end
    %load already generated tempfile
    load(tempfile,'outputData','isNoise');
    [~,tempbase,~] = fileparts(tempfile);
    tempsplit = strsplit(tempbase,'_');
    tempsuff = strjoin(['temp' tempsplit(1:(length(tempsplit)-1))],'_');
    valsuff = [varpref num2str(varval)];
    
    %also load options!
    load(tempfile,'options');
    options.(varname) = varval;

    %Now assign the validation data to templates
    [outdir2,outbase2,~] = fileparts(valsignaldir);
    outdir2 = [outdir2 '/output_' outbase2];
    outbase2 = strjoin({outbase2, tempsuff, valsuff},'_');
    if exist(outdir2,'dir') ~= 7
        mkdir(outdir2);
    end
    %Loop over validation signal files
    signallist = dir([valsignaldir '/*.wav']);
    noiselist = dir([valnoisedir '/*.wav']);
    signalnames = cell(1,length(signallist));
    noisenames = cell(1,length(noiselist));
    for s=1:length(signallist)
        signalnames{s}=signallist(s).name;
    end
    for n=1:length(noiselist)
        noisenames{n}=noiselist(n).name;
    end
    signalnames = signalnames(~cellfun('isempty',signalnames));
    noisenames = noisenames(~cellfun('isempty',noisenames));
    signalThresholds = [];
    signalAssignments = {}; %save filenames to this cell array
    %Instead, create a new folder within the original, and save the
    %clustering output here
    outdirsignal = fullfile(valsignaldir,['clusteringoutput_' tempsuff]);
    if exist(outdirsignal,'dir') ~= 7
        mkdir(outdirsignal);
    end
    for s=1:length(signalnames)
        isBadFile = strfind(signalnames{s},'._');
        if isempty(isBadFile)
            [dfilt2,fs2] = audioread(fullfile(valsignaldir,signalnames{s}));
            options.fs = fs2;
            [groupings,peakIdxGroup,likes,allPeakIdx,allNormalizedPeaks,noiseThreshold] = ...
                     assignDataToTemplates(dfilt2,outputData,options);
            %Add outputs to a structure? This seems very tricky to save.
            [~,sname,~]=fileparts(signalnames{s});
            outbase2 = strjoin({sname, tempsuff, valsuff},'_');
            save(fullfile(outdirsignal,[outbase2 '_outputAssignTemplates.mat']), 'groupings',...
        'peakIdxGroup','likes','allPeakIdx','allNormalizedPeaks','noiseThreshold');
            %outAssign = struct('groupings',groupings,...
             %   'peakIdxGroup',peakIdxGroup,'likes',likes,'allPeakIdx',...
              %  allPeakIdx,'allNormalizedPeaks',allNormalizedPeaks,...
               % 'noiseThreshold',noiseThreshold);
            %outAssign.peakIdxGroup = peakIdxGroup;
            signalAssignments{s} = fullfile(outdirsignal,[outbase2 '_outputAssignTemplates.mat']);
            signalThresholds = [signalThresholds noiseThreshold];
        end
    end
    %Don't make any plots
    
    %Now run for noise files
    %fprintf(['The minimum noiseThreshold for signal files (to be used for\n' ...
    %    'min_noise_threshold) is %.3f (converts to %.3f) for %s of %.2f\n'],...
    %    min(signalThresholds),sqrt(10^min(signalThresholds)),varname,varval);
    options.min_noise_threshold = sqrt(10^min(signalThresholds));
    noiseAssignments = {};
    outdirnoise = fullfile(valnoisedir,['clusteringoutput_' tempsuff]);
    if exist(outdirnoise,'dir') ~= 7
        mkdir(outdirnoise);
    end
    for n=1:length(noisenames)
        isBadFile = strfind(noisenames{n},'._');
        if isempty(isBadFile)
            [dfilt2,fs2] = audioread(fullfile(valnoisedir,noisenames{n}));
            options.fs = fs2;
            [groupings,peakIdxGroup,likes,allPeakIdx,allNormalizedPeaks,noiseThreshold] = ...
                     assignDataToTemplates(dfilt2,outputData,options);
            %Add outputs to a structure?
            [~,nname,~]=fileparts(noisenames{n});
            outbase2 = strjoin({nname, tempsuff, valsuff},'_');
            save(fullfile(outdirnoise,[outbase2 '_outputAssignTemplates.mat']), 'groupings',...
        'peakIdxGroup','likes','allPeakIdx','allNormalizedPeaks','noiseThreshold');
            %outAssign = struct('groupings',groupings,...
             %   'peakIdxGroup',peakIdxGroup,'likes',likes,'allPeakIdx',...
              %  allPeakIdx,'allNormalizedPeaks',allNormalizedPeaks,...
               % 'noiseThreshold',noiseThreshold);
            %outAssign.peakIdxGroup = peakIdxGroup;
            noiseAssignments{n} = fullfile(outdirnoise,[outbase2 '_outputAssignTemplates.mat']);
        end
    end
    %Save each set of outputs together
    %Redo the name making
    [outdir2,outbase2,~] = fileparts(valsignaldir);
    outdir2 = [outdir2 '/output_' outbase2];
    outbase2 = strjoin({outbase2, tempsuff, valsuff},'_');
    save(fullfile(outdir2,[outbase2 '_outputAssignTemplates.mat']), ...
        'signalAssignments','noiseAssignments','signalnames',...
        'noisenames','options');
    
    %Modify EstimateFPFN to work with signalAssignments and
    %noiseAssignments
    [~,~,~] = EstimateFPFN_SeparateSongs(tempfile,fullfile(outdir2,...
        [outbase2 '_outputAssignTemplates.mat']),categorytimes);
    fprintf(['The minimum noiseThreshold for signal files (to be used for\n' ...
        'min_noise_threshold) is %.3f (converts to %.3f) for %s of %.2f\n'],...
        min(signalThresholds),sqrt(10^min(signalThresholds)),varname,varval);
    %Write filename from EstimateFPFN to filelist
    fpfnfile = fullfile(outdir2,[outbase2 '_outputAssignTemplates_fpfn.mat']);
    flID = fopen(filelist,'a');
    fprintf(flID,'%s\n',fpfnfile);
    fclose(flID);
end