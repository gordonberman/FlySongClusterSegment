function randomizeAndCombineSignalNoiseWav(signalfolder, noisefolder, outname, outfs)
    %List all .wav files in signalfolder and noisefolder
    %Combine them in a random order
    %Write as a single .wav file (outname.wav)
    %with sampling frequency outfs
    %Write file order (outname.order.csv)'
    signallist = dir([signalfolder '/*.wav']);
    noiselist = dir([noisefolder '/*.wav']);
    signalnames = cell(1,length(signallist));
    noisenames = cell(1,length(noiselist));
    for s=1:length(signallist)
        signalnames{s}=signallist(s).name;
    end
    for n=1:length(noiselist)
        noisenames{n}=noiselist(n).name;
    end
    signalnames(~cellfun('isempty',signalnames));
    noisenames(~cellfun('isempty',noisenames));
    allnames=horzcat(signalnames,noisenames);
    neworder=randperm(length(allnames));
    wavlist=cell(1,length(allnames));
    for w=1:length(allnames)
        wavlist{w} = allnames{neworder(w)};
    end
    %write wavlist to a file
    %if possible, also add signal/noise, as well as length
    %and cumulative length
    %formatSpec = '%s\n';
    %fileID = fopen([outname '.txt'],'w');
    %for f = 1:length(wavlist)
    %    fprintf(fileID,formatSpec,wavlist{f});
    %end
    %fclose(fileID);
    category=cell(1,length(allnames));
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
            %determine whether file belongs to signalfolder or noisefolder
            if (sum(strcmp(signalnames,wavlist{w}))>0)
                wavfolder=signalfolder;
                category{w}='signal';
            elseif (sum(strcmp(noisenames,wavlist{w}))>0)
                wavfolder=noisefolder;
                category{w}='noise';
            else
                error('Issue with file: %s\n',wavlist{w});
            end
            [inwav, infs] = audioread([wavfolder '/' wavlist{w}]);
            if infs ~= outfs
                sprintf(['Sampling frequency for %s (%d) not equal to '...
                         'desired sampling frequency (%d); converting.'],...
                         wavlist{w}, infs, outfs);
                [P,Q] = rat(outfs/infs);
                inwav = resample(inwav,P,Q);
            end
            timecounter = timecounter + length(inwav); %fixed so length is after resampling
            cumtime{w} = timecounter;
            outwav = vertcat(outwav,inwav);
        end
    end
    wavlist=wavlist';
    category=category';
    cumtime=cumtime';
    outtable = table(wavlist,category,cumtime);
    writetable(outtable, [outname '.order.csv']);
    audiowrite([outname '.wav'], outwav, outfs);
end