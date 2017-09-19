function wrapCreateTraining(folderorlist,...
    goalLengthInSPerSpecies,repopath,minFromEachDataSet)
    %Run createTrainingSet for set of wav files in folder, list, or listoflists
    %Randomize and combine the resulting outputs
    %Input:
    %folderorlist (character): one of the following:
    %1) a folder containing .wav files to be used in the training set
    %2) a file listing the .wav files to be used in the training set
    %3) a file listing lists as in (2), one for each species
    %Output (files written):
    %fullfile(fldir,[flbase '_trainingset.wav']): training set (wav file)
    %fullfile(fldir,[flbase '_trainingset.ordered.csv']): where the
    %training set came from
    %fullfile(fldir,[flbase '_thresholds.csv']): per-file noise thresholds

    %Add paths:
    %Set some defaults:
    if nargin < 3 || isempty(goalLengthInSPerSpecies)
        goalLengthInSPerSpecies = 300; %How many seconds per species?
    end
    if nargin < 4 || isempty(repopath)
        repopath = '/net/mahler/wynn/repos/FlySongClusterSegment';
    end
    if nargin < 5 || isempty(minFromEachDataSet)
        minFromEachDataSet = 0; %change this if all of your files definitely 
                                %have song
    end
    addpath(genpath(repopath));
    
    listofwav = true; %change to false if it's a list of lists
    if exist(folderorlist,'dir') == 7
        wavfolder = folderorlist;
        %determine outfile name
        [fldir,flbase,~]=fileparts(wavfolder);
        %get list of all  .wav files in the folder
        dirwavlist = dir([wavfolder '/*.wav']); %will contain some '.' files
        for k = length(dirwavlist):-1:1
            fname = dirwavlist(k).name;
            if fname(1) == '.'
                dirwavlist(k) = [];
            end
        end
        allfiles = {length(dirwavlist)};
        for f = length(dirwavlist):-1:1
            allfiles{f} = fullfile(dirwavlist(f).folder,dirwavlist(f).name);
        end
    elseif exist(folderorlist,'file') == 2
        wavlist = folderorlist;
        %read in file
        [fldir,flbase,~]=fileparts(wavlist);
        folfileID=fopen(wavlist);
        folders=textscan(folfileID,'%s');
        fclose(folfileID);
        allfiles=folders{1}; %scans as a cell array of cell arrays
        [~,~,firstext] = fileparts(allfiles{1});
        if strcmp(firstext,'.txt')
            listofwav = false;
        elseif strcmp(firstext,'.wav') == 0
            error('File in %s is not a recognized type\n', allfiles{1});
        end
    else
        error('%s does not exist!\n',folderorlist);
    end
    
    %loop over listoflists
    if listofwav
        allfiles = {allfiles};
    end
    convertrates = false;
    fs = cell(length(allfiles),1);
    speciesfiles = cell(length(allfiles),1);
    containsmat = false;
    for a=1:length(allfiles)
        %check for equivalent fs
        if ~listofwav
            afID=fopen(allfiles{a});
            speciesfilescell=textscan(afID,'%s');
            fclose(afID);
        else
            speciesfilescell=allfiles;
        end
        speciesfiles{a}=speciesfilescell{1};
        if exist('catspeciesfiles','var') ~= 1
            catspeciesfiles = speciesfiles{a};
        else
            catspeciesfiles = [catspeciesfiles; speciesfiles{a}];
        end
        fs{a} = cell(length(speciesfiles{a}),1);
        for s=1:length(speciesfiles{a})
            %Allow files to be .mat files with data.fs field
            [~,fname,fext] = fileparts(speciesfiles{a}{s});
            if strcmp(fext,'.wav')
                try
                    sinfo=audioinfo(speciesfiles{a}{s});
                    fs{a}{s} = sinfo.SampleRate;
                catch MEnotfound
                    speciesfiles{a}{s}
                    rethrow(MEnotfound);
                end
            elseif strcmp(fext,'.mat')
                containsmat = true;
                try
                    sdata = load(speciesfiles{a}{s},'data');
                    fs{a}{s} = sdata.data.fs;
                catch MEfs
                    error('No data.fs found in %s\n',[fname fext]);
                end
            else
                error('File %s is not a recognized type\n', [fname fext]);
            end
            if exist('minfs','var') ~= 1
                minfs = fs{a}{s};
            else
                if fs{a}{s} ~= minfs
                    convertrates = true;
                    minfs = min(minfs,fs{a}{s});
                end
            end
        end
    end
    
    %get trainingSetLength using minfs
    trainingSetLength = goalLengthInSPerSpecies*minfs;
    
    %save combined output, and scramble order afterwards
    for f=1:length(allfiles)
        %fprintf(1, 'Creating training set for %s\n',allfiles{f});
        fprintf(1, 'Creating training set for set %i of %s\n',f,folderorlist);
        %run createTrainingSet
        if convertrates || containsmat
            fprintf(1,['   %s either contains some files with\n' ...
                '   unequal sampling rate or contains non-wav files\n'],...
                flbase);
            data = {};
            for l=1:length(speciesfiles{f})
                [~,fname,fext] = fileparts(speciesfiles{f}{l});
                if strcmp(fext,'.wav')
                    [tmpsong,tmpfs]=audioread(speciesfiles{f}{l});
                elseif strcmp(fext,'.mat')
                    try
                        sdata = load(speciesfiles{f}{l},'data');
                        tmpsong = sdata.data.d;
                        tmpfs = sdata.data.fs;
                    catch MEfs
                        error('No data.fs or data.d found in %s\n',[fname fext]);
                    end
                else
                    error('File %s is not a recognized type\n', [fname fext]);
                end
                if tmpfs ~= minfs
                    tmpsong = fixSamplingFrequency(tmpsong,tmpfs,minfs);
                end
                data = [data; tmpsong];
            end
        else
            data = speciesfiles{f};
        end
        options.fs = minfs; %also possible to modify other options
        [mintrainingSet,mintrainingSet_origins,minthresholds] = ...
            createTrainingSet_wm111116(data,trainingSetLength,...
            minFromEachDataSet,options);
        fprintf(1,'%.0f snippets included from set %i of %s\n',...
            size(mintrainingSet_origins,1),f,folderorlist);
        %fprintf(1,'%.0f snippets included from %s\n',...
        %    size(mintrainingSet_origins,1),allfiles{f});
        %Convert to table and add filename BEFORE combining
        mininfotab = array2table(mintrainingSet_origins,...
            'VariableNames',{'Start','End','Dataset','RecStart','RecEnd'});
        mininfotab.File = num2cell(mininfotab.Dataset);
        for m=1:height(mininfotab)
            mininfotab.File{m} = speciesfiles{f}{mininfotab.Dataset(m)};
        end
        if exist('basetrainingSet','var') ~= 1
            basetrainingSet = mintrainingSet;
            %basetrainingSet_origins = mintrainingSet_origins;
            baseinfotab = mininfotab;
            thresholds = minthresholds;
        else
            basetrainingSet = [basetrainingSet; mintrainingSet];
            %Fix Start and End
            lastend = max(baseinfotab.End);
            mininfotab.Start = mininfotab.Start + lastend;
            mininfotab.End = mininfotab.End + lastend;
            baseinfotab = [baseinfotab; mininfotab];
            thresholds = [thresholds; minthresholds];
        end
    end
    
    %combine outputs, randomizing order
    fprintf(1,'   Randomizing order for training set from %s\n',flbase);
    trainingSet = zeros(length(basetrainingSet),1);
    traininfotab = baseinfotab;

    tord = randperm(height(baseinfotab));
    starttset = 1;
    for t=1:length(tord)
        traininfotab(t,:) = baseinfotab(tord(t),:);
        lengthtoadd = traininfotab.End(t) - traininfotab.Start(t);
        trainingSet(starttset:(starttset+lengthtoadd)) = ...
            basetrainingSet(traininfotab.Start(t):traininfotab.End(t));
        traininfotab.Start(t) = starttset;
        traininfotab.End(t) = starttset+lengthtoadd;
        starttset = starttset + lengthtoadd + 1;
    end

    %write output files
    %[fldir,flbase,~]=fileparts(listoflists);
    audiowrite(fullfile(fldir,[flbase '_trainingset.wav']),trainingSet,minfs);
    writetable(traininfotab,fullfile(fldir,[flbase '_trainingset.ordered.csv']));
    %also write noise thresholds
    convthresholds = sqrt(10.^thresholds);
    fprintf(1,['The minimum noise threshold for files in the training ' ...
        'set is %.3f, or %.3f when converted; the latter should be ' ...
        'used to set minnoise.'], min(thresholds), min(convthresholds));
    if size(catspeciesfiles,2) > size(catspeciesfiles,1)
        catspeciesfiles = catspeciesfiles';
    end
    fprintf(1,'Creating threshold table\n');
    threshtab = table(catspeciesfiles,thresholds,convthresholds);
    fprintf(1,'Writing thresholds\n');
    writetable(threshtab,fullfile(fldir,[flbase '_thresholds.csv']));
    fprintf(1,'Finished creating training set\n');
end
