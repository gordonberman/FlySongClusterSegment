function randomizeWavForTraining(wavfolder, outname, outfs)
    %List all .wav files in wavfolder
    %Combine them in a random order
    %Write as a single .wav file (outname.wav)
    %with sampling frequency outfs
    %Write file order (outname.order.csv)
    if ~isnumeric(outfs)
        outfs = str2num(outfs);
    end
    signallist = dir([wavfolder '/*.wav']);
    allnames = cell(1,length(signallist));
    for s=1:length(signallist)
        allnames{s}=signallist(s).name;
    end
    allnames(~cellfun('isempty',allnames));
    neworder=randperm(length(allnames));
    wavlist=cell(1,length(allnames));
    for w=1:length(allnames)
        wavlist{w} = allnames{neworder(w)};
    end
    if length(wavlist) == 0
        error('No wav files found in %s!\n',wavfolder);
    end
    %write wavlist to a file
    %also add cumulative length

    cumtime=cell(1,length(allnames));
    timecounter=0;
    %concatenate wav files from wavlist
    for w=1:length(wavlist)
        if exist('outwav','var')~=1
            outwav=[];
            category{w}='NA';
            cumtime{w}='NA';
        end
        isBadFile = strfind(wavlist{w},'._');
        if isempty(isBadFile)
            [inwav, infs] = audioread([wavfolder '/' wavlist{w}]);
            timecounter = timecounter + length(inwav);
            cumtime{w} = timecounter;
            if infs ~= outfs
                sprintf(['Sampling frequency for %s (%d) not equal to '...
                         'desired sampling frequency (%d); converting.'],...
                         signallist(w).name, infs, outfs);
                [P,Q] = rat(outfs/infs);
                inwav = resample(inwav,P,Q);
            end
            outwav = vertcat(outwav,inwav);
        end
    end
    wavlist=wavlist';
    cumtime=cumtime';
    outtable = table(wavlist,cumtime);
    writetable(outtable, [outname '.order.csv']);
    audiowrite([outname '.wav'], outwav, outfs);
end
