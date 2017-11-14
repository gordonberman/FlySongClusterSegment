function multiSpeciesCreateTraining(listoflists)
    %Run createTrainingSet for each filelist in listoflists
    %Randomize and combine the resulting outputs
    %Input:
    %listoflists (character): a file listing each species' recording list
    %Output (files written):
    %Add paths:
    addpath(genpath('/global/scratch/wynn/repos/bachtroglabsong'),...
            genpath('/global/scratch/wynn/repos/FlySongClusterSegment'));
    
    %Set some defaults:
    goalLengthInSPerSpecies = 60; %How many seconds per species?
    minFromEachDataSet = 0; %change this if all of your files definitely 
                            %have song
    
    if exist(listoflists,'file') ~= 2
        error('%s does not exist!\n',listoflists);
    end
    %read in listoflists
    [fldir,flbase,~]=fileparts(listoflists);
    folfileID=fopen(listoflists);
    folders=textscan(folfileID,'%s');
    fclose(folfileID);
    allfiles=folders{1}; %scans as a cell array of cell arrays
    
    %loop over listoflists
    convertrates = false;
    fs = cell(length(allfiles),1);
    speciesfiles = cell(length(allfiles),1);
    for a=1:length(allfiles)
        %check for equivalent fs
        afID=fopen(allfiles{a});
        speciesfilescell=textscan(afID,'%s');
        fclose(afID);
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
                sinfo=audioinfo(speciesfiles{a}{s});
                fs{a}{s} = sinfo.SampleRate;
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
        fprintf(1, 'Creating training set for %s\n',allfiles{f});
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
        fprintf(1,'%.0f snippets included from %s\n',...
            size(mintrainingSet_origins,1),allfiles{f});
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
    threshtab = table(catspeciesfiles,thresholds,convthresholds);
    writetable(threshtab,fullfile(fldir,[flbase '_thresholds.csv']));
end
