function GetSongNoiseLevels(songlist,maxmaxP,fscsdir)
    %songlist: text file listing full path to song files to evaluate
    %maxmaxP: maxP cutoff for files to go in clean list
    %fscsdir: directory for FlySongClusterSegment
    %to do: get 'threshold' as in createTrainingSet
    %add paths
    addpath(genpath(fscsdir));
    %generate options (can use an options file from createTemplates if available)
    options.fs=6000;
    if ischar(maxmaxP)
        maxmaxP = str2num(maxmaxP);
    end
    %load file list
    fullfid=fopen(songlist,'r');
    fullfilescell = textscan(fullfid,'%s','Delimiter','\n');
    fullfiles = fullfilescell{1};
    fclose(fullfid);
    filestosave = fullfiles; %exclude ones with high maxP
    %loop over files, getting maxP for each
    maxP=zeros(length(fullfiles),1);
    sdNoise=zeros(length(fullfiles),1);
    threshold=zeros(length(fullfiles),1);
    for f=1:length(fullfiles)
        if exist(fullfiles{f},'file')==2
            [~,~,fext] = fileparts(fullfiles{f});
            if strcmp(fext,'.wav')
                [song,Fs]=audioread(fullfiles{f});
            elseif strcmp(fext,'.mat')
                data = load(fullfiles{f});
                song = data.data.d;
                Fs = data.data.fs;
            else
                error('File %s is not in a recognizable format\n',fullfiles{f});
            end
            [tmpnoiseData,~,~,maxP(f),~,threshold(f),~]=createNoiseDataSet(song,Fs*300,options);
            sdNoise(f) = std(tmpnoiseData);
            if maxP(f) > maxmaxP
                filestosave{f} = 'na'; %somehow eliminate
            end
        else
            fprintf(1,'  Cant find %s\n',fullfiles{f});
            continue
        end
    end
    filestosave(find(strcmp(filestosave,'na'))) = [];
    %make a table with filenames and maxP values
    convthresholds = sqrt(10.^threshold);
    fprintf(1,['The minimum noise threshold for files in the training ' ...
        'set is %.3f, or %.3f when converted; the latter should be ' ...
        'used to set minnoise.'], min(threshold), min(convthresholds));
    outtab = table(fullfiles,maxP,sdNoise,threshold,convthresholds);
    [filedir,filebase,~]=fileparts(songlist);
    writetable(outtab,fullfile(filedir,[filebase '_withmaxP.csv']));
    %also write list of files under maxmaxP threshold
    fileID = fopen(fullfile(filedir,[filebase '_cleansongs.txt']),'w');
    for f=1:length(filestosave)
    	fprintf(fileID,'%s\n',filestosave{f});
    end
    fclose(fileID);
end
