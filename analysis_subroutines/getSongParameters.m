function [pulsetrains, options] = getSongParameters(data, peakIdxGroup, ...
    isNoise, dogrouped, options, freqIdxGroup)
    %processing script downstream of assignDataToTemplates
    %pulsetrains should be a structure containing data
    %on each pulse train (location of all pulses and their CF,
    %all IPI, PTL, PPB, templateGroups)
    %necessary options: fs, distance separating pulses
    %to be considered a new pulse train
    %reporting pulses, ipi, and ptl in milliseconds, cf in Hz
    
    %Get the full options (including numIPIBins and IPI_sigma)
    %The two above are set in oneSongSummarize....
    options.setAll = false;
    options = makeParameterStructure(options);
    %If times are in ms, don't convert IPI_sigma
    %IPI_sigma = options.IPI_sigma * options.fs / 1000;
    
    %only do the following if freqIdxGroup not provided:
    signalPeakIdx = getSignalPeakIdx(peakIdxGroup, isNoise);
    %need to find carrier frequency just for signal peaks as well
    %if carrier frequency specified:
    %[signalPeakIdx, signalFreqIdx] = getSignalFreqIdx(freqIdxGroup, isNoise);
    if dogrouped
        signalPeakGroup = getSignalPeakGroup(signalPeakIdx, peakIdxGroup, ...
            isNoise);
        nagroup = length(peakIdxGroup) + 1;
    else
        signalPeakIdx = getSignalPeakIdx(peakIdxGroup, isNoise);
        signalPeakGroup = zeros(length(signalPeakIdx));
        nagroup = 1;
    end
    ts = signalPeakIdx*1000/options.fs; %convert to milliseconds
    ipiall = diff(ts);
    
    %find max of IPI distribution through kernel density estimation
    [Y,X] = hist(ipiall,options.numIPIBins);
    Y = normalizeHist(X,Y);
    IPI_sigma = options.IPI_sigma / (X(2) - X(1)); %convert to units of bin widths
    Y = gaussianfilterdata(Y,IPI_sigma);   
    s = fit(X',Y','spline');
    modeIPI_estimate = fminsearch(@(x) -s(x),X(argmax(Y)));
    %Use the maximum modeIPI_estimate across recordings to set summmaxIPI:
    %summmaxIPI = modeIPI_estimate*3

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
        tmpgroup = [signalPeakGroup(CC.PixelIdxList{i}); ...
            signalPeakGroup(CC.PixelIdxList{i}(end)+1)];
        [groupmode,~,withmult] = mode(tmpgroup);
        if length(withmult) > 1
            templateGroups(i) = nagroup;
        else
            templateGroups(i) = groupmode;
        end
        ipi{i} = diff(pulses{i});
        %get carrier frequency if available
        %carrier frequency estimation will require some distance around
        %peak?
        %cf{i} = getCarrierFrequency(pulses{i});
    end
    pulsetrains = struct('pulses',pulses,'ipi',ipi,'cf',cf,'numPulses',...
        numPulses,'trainLengths',trainLengths,'templateGroups',...
        templateGroups);
end