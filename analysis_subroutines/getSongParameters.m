function [pulsetrains, options] = getSongParameters(wavfile, peakIdxGroup, ...
    isNoise, dogrouped, options, freqIdxGroup)
    %processing script downstream of assignDataToTemplates
    %pulsetrains should be a structure containing data
    %on each pulse train (location of all pulses and their CF,
    %all IPI, PTL, PPB, templateGroups)
    %necessary options: fs, distance separating pulses
    %to be considered a new pulse train
    %reporting pulses, ipi, and ptl in milliseconds, cf in Hz
    %first option used to be data, now wavfile:
    if isnumeric(wavfile)
        fprintf(1,'Getting pulse trains for data of length %d\n',...
            length(wavfile));
    else
        fprintf(1,'Getting pulse trains for %s\n',wavfile);
    end
    
    %Get the full options (including numIPIBins and IPI_sigma)
    %The two above are set in oneSongSummarize....
    options.setAll = false;
    options = makeParameterStructure(options);
    %If times are in ms, don't convert IPI_sigma
    %IPI_sigma = options.IPI_sigma * options.fs / 1000;
    
    %only do the following if freqIdxGroup not provided:
    if nargin < 6 || isempty(freqIdxGroup)
        signalPeakIdx = getSignalPeakIdx(peakIdxGroup, isNoise);
    else
        %need to find carrier frequency just for signal peaks as well
        %if carrier frequency specified:
        [signalPeakIdx, signalFreqIdx] = getSignalFreqIdx(peakIdxGroup,...
            freqIdxGroup, isNoise);
    end
    if dogrouped
        signalPeakGroup = getSignalPeakGroup(signalPeakIdx, peakIdxGroup, ...
            isNoise); %groups pertain to both peak Idx and Freq
        nagroup = length(peakIdxGroup) + 1;
    else
        %signalPeakIdx = getSignalPeakIdx(peakIdxGroup, isNoise);
        signalPeakGroup = zeros(1,length(signalPeakIdx));
        nagroup = 1;
    end
    ts = signalPeakIdx*1000/options.fs; %convert to milliseconds
    ipiall = diff(ts);

    %Now use summmaxIPI to find pulse trains
    isConnected = ipiall <= options.summmaxIPI;
    CC = bwconncomp(isConnected);
    pulses = cell(1, CC.NumObjects);
    cf = cell(1, CC.NumObjects);
    ipi = cell(1, CC.NumObjects);
    numPulses = cell(1, CC.NumObjects);
    trainLengths = cell(1, CC.NumObjects);
    templateGroups = zeros(1, CC.NumObjects);
    for i=1:CC.NumObjects
        numPulses{i} = length(CC.PixelIdxList{i});
        trainLengths{i} = ts(CC.PixelIdxList{i}(end)+1) - ts(CC.PixelIdxList{i}(1));
        %pulses are now in milliseconds
        pulses{i} = [ts(CC.PixelIdxList{i}) ts(CC.PixelIdxList{i}(end)+1)];
        %the same indexing should also work for carrier frequencies
        if exist('signalFreqIdx','var') == 1
                cf{i} = [signalFreqIdx(CC.PixelIdxList{i}) ...
                    signalFreqIdx(CC.PixelIdxList{i}(end)+1)];
        end
        tmpgroup = [signalPeakGroup(CC.PixelIdxList{i}) ...
            signalPeakGroup(CC.PixelIdxList{i}(end)+1)];
        [groupmode,~,withmult] = mode(tmpgroup);
        if length(withmult) > 1
            templateGroups(i) = nagroup;
        else
            templateGroups(i) = groupmode;
        end
        ipi{i} = diff(pulses{i});
    end
    pulsetrains = struct('pulses',pulses,'ipi',ipi,'cf',cf,'numPulses',...
        numPulses,'trainLengths',trainLengths,'templateGroups',...
        templateGroups);
end
